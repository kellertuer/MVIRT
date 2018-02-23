function eta = gradData(M,f,x,p,epsilon)
% gradData(M,f,x)
% gradient for 1/p d^p(x,f)
if nargin < 4
    p=2;
end
if nargin < 5
    epsilon=0;
dataDim = size(f);
dataDim = dataDim(  (length(M.ItemSize)+1) : end  );
if p ~= 2
    d = M.dist(f,x);
    eta = -M.log(x,f);
    eta = eta.*permute(d.^(p-2),[(length(dataDim)+1):length(size(f)), 1:length(dataDim)]);
else
    eta = -M.log(x,f);
end
end

