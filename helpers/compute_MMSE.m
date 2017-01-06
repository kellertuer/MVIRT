function [u_estimate,mu,Sigma] = compute_MMSE(M,u,u_oracle,sigma,step,simple)
% [u_estimate,mu,Sigma] = compute_MMSE(M,u,u_oracle,sigma,step,simple)
% compute an estimate based on a second order approach wrt to the original
% and eventually oracle image
%
% INPUT:
%   M           Manifold
%   u           noisy image patches  manDim x patch_size x patch_size x K
%   u_oracle    denoised image patches (size as u)
%   sigma       standard deviation of noise
%   step        indicates which step of the algorihtm we are in
%   simple      indicates if only means are used
% OUTPUT:
%   u_estimate  restored patches
%   mu          mean patch
%   Sigma       emperical covariance matric of the patch
%
%   For further details see
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%     Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
% ---
% Manifold-valued Image Restoration Toolbox 1.2
%  J. Persch  ~2017-07-05 | R. Bergmann 2016-01-07
% see LICENSE.txt

dimen = size(u);
patch_size = dimen(end-2);
manDim = M.ItemSize;
coef = M.Dimension;
K = dimen(end);
lIS = length(manDim);
L = patch_size^2;
% Reshape patches
observations = reshape(u,[M.ItemSize,L,K]);
mu = M.mean(observations);
if ~simple % Second Order Update
    observations_ref = reshape(u_oracle,[M.ItemSize,L,K]);
    % Shift Observation to Tangential Plane
    observation = M.log(repmat(mu,[ones(1,lIS+1),K]),observations);
    observation_ref = M.log(repmat(mu,[ones(1,lIS+1),K]),observations_ref);
    
    
    % Find basis representation for observation
    ONB = M.TpMONB(mu);
    % Calculate the coefficients of the observations
    coef_obs = permute(M.dot(repmat(mu,[ones(1,lIS+1),K,coef]),repmat(permute(ONB,[1:lIS+1  lIS+3 lIS+2 ]),[ones(1,lIS+1),K,1]),...
        repmat(observation,[ones(1,length(dimen)-1),coef])),[3,1,2]);
    coef_obs = reshape(coef_obs,[],K);
    
    coef_obs_ref = permute(M.dot(repmat(mu,[ones(1,lIS+1),K,coef]),repmat(permute(ONB,[1:lIS+1  lIS+3 lIS+2 ]),[ones(1,lIS+1),K,1]),...
        repmat(observation_ref,[ones(1,length(dimen)-1),coef])),[3,1,2]);
    coef_obs_ref = reshape(coef_obs_ref,[],K);
    
    Sigma = zeros(coef*L);
    % Calculate covariance matrix
    for i = 1:K
        Sigma = Sigma+coef_obs_ref(:,i)*coef_obs_ref(:,i)';
    end
    Sigma = 1/(K-1)*Sigma;
    % Perform covariance update
    if cond(Sigma) < 10^6
        if step == 1
            A = (Sigma-sigma^2*eye(coef*L))/Sigma;
        else
            A = Sigma/(Sigma+sigma^2*eye(coef*L));
        end
        coef_obs = A*coef_obs;
    else
        Sigma = zeros(size(Sigma));
    end
    % Put the new coefficients back onto the manifold and reshape
    ONB = permute(reshape(ONB,prod(manDim),L,coef),[1,3,2]);
    u_estimate = zeros(prod(manDim),K,L);
    for i = 1:L
        u_estimate(:,:,i) = ONB(:,:,i)*coef_obs((i-1)*coef+(1:coef),:);
    end
    u_estimate = reshape(M.exp(repmat(mu,[ones(1,lIS+1),K]),reshape(permute(u_estimate,[1,3,2]),[manDim,L,K])),dimen);
else
    u_estimate = repmat(mu,[ones(1,lIS+1),K]);
    Sigma = 1;
end
mu = reshape(mu,dimen(1:end-1));
end

