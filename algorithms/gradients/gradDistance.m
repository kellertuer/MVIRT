function eta = gradDistance(M,x,f,p)
% gradData(M,f,x)
% gradient for 1/p d^p(x,f)
if nargin < 4
    p=2;
end
mD = length(M.ItemSize);
if p ~= 2
    d = M.dist(f,x);
    eta = -M.log(x,f);
    eta = eta.*shiftdim(d.^(p-2),-mD);
    % securly avoid division by zero
    if any(d(:)<eps)
        warning('gradDist returns an element from the subdifferential only.');
        eta(M.manDims{:},d<eps) = 0;
    end
else
    eta = -M.log(x,f);
end
end