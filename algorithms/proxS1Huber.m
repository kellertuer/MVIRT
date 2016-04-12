function x = proxS1Huber(pf,lambda,tau, omega)
%proxacdHuber(f, lambda,q,w) compute the [PROX]imal mapping for
%       [A]bsolute [C]yclic [D]ifferences of first order, where the Huber
%       functional is employed using tau and omega, i.e. for one data set
%       f_1,f_2 we employ h(d(f_1,f_2)), where h(s) = 
%       tau^s^ for s < omega/(sqrt(2)*tau)
%       omega*sqrt(2)tau s - \omega^/2 else. 
%
%   INPUT
%       f      : each column of f is threatened as one data set
%                  if f is a row vector, it isjust threated as one data set
%       lambda : usual value of the proximum
%       tau    : 
%       omega  : 
%
%   OUTPUT
%       x      : proximal mappings
%
% SEE ALSO
% 
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014). 
%       Second order differences of cyclic data and applications in variational denoising. 
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
%    and
%    A. Weinmann, L. Demaret, M. Storath:
%       Total Variation Regularization for Manifold-valued Data
%       SIAM Journal on Imaging Sciences. 7, (4), 2226?2257.
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-08-27 | 2016-02-22
% see LICENSE.txt

if isrow(pf)
    f = pf';
else
    f = pf;
end
[n,m] = size(f); %number of sets (rows, m) and length of each data set (columns, n)
if (n~=2)
    error('The length of each data set has to be 2.');
end
if (omega <=0) || (tau <=0)
    error('This method only provides proximal mappings omega, tau > 0.');
end
% and the rest is not handled at all
sp = symMod( dot( repmat([-1,1]',[1,m]),f), 2*pi);
sc = abs(sp) < omega*(1+4*lambda*tau^2)/(sqrt(2)*tau); %small cases, i.e. squared
% By Theorem 3.5 of BSW14 - ommitting the + case for pi
s = sign(sp);
% elementwise minima, denoted m in Thm 3.5 and adapted to huber
mins = min(ones(1,m).*lambda*sqrt(2)*omega*tau, abs(sp)/2); %outer case similar to q=1
% Special ?small? cases - employ huber quadratic case
mins(sc) = 2*lambda*tau^2/(1+4*lambda*tau^2);
% Multiply elementwise the sign (s) - the rest works the same as
% usual mult
x = symMod(f - [-1,1]'*(s.*mins),2*pi); % f + s*m*w for each row of f
if isrow(pf)
    x = x';
end
end

