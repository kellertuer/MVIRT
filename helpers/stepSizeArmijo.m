function stepsize = stepSizeArmijo(varargin)

ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'F');
addRequired(ip,'x');
addRequired(ip,'descentDir');
addOptional(ip,'Gradient',[]);% if not set assume descentDir = -grad
addOptional(ip,'InitialStepsize',1);
addOptional(ip,'rho',0.5);
addOptional(ip,'c',0.0001);
parse(ip, varargin{:});
vars = ip.Results;
if isempty(vars.Gradient)
    Gradient = -vars.descentDir;
elseif all(size(vars.Gradient) == size(vars.eta))
    Gradient = vars.Gradient;
else
    error('Gradient and descent direction have to be of same dimension');
end

e=vars.F(vars.x);
stepsize = vars.InitialStepsize;
e_n = e-1;
while e > e_n
    x_n = vars.M.exp(vars.x,vars.descentDir,stepsize);
    e_n = vars.F(x_n)+vars.c*stepsize*sum(...
        vars.M.dot(vars.x(vars.M.allDims{:},:),-vars.descentDir(vars.M.allDims{:},:),Gradient(vars.M.allDims{:},:)));
    stepsize = stepsize/vars.rho;
end
stepsize = stepsize*vars.rho;
while e < e_n && stepsize > eps
    stepsize = stepsize*vars.rho;
    x_n = vars.M.exp(vars.x,vars.descentDir,stepsize);
    e_n = vars.F(x_n)+vars.c*stepsize*sum(...
        vars.M.dot(vars.x(vars.M.allDims{:},:),-vars.descentDir(vars.M.allDims{:},:),Gradient(vars.M.allDims{:},:)));
end
