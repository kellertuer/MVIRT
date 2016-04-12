function s = TeXNum(x,fmtmant,fmtexpn)
% TeXNum(x,fmt1,fmt2) 
%    produce a number in scientific notation with respect to base 10, where
%
%  INPUT
%    x       : a number
%    fmtmant : a number format for the mantisse, e.g. '%1.2f'
%    fmtexpn : a number format for the exponent, e.g. '%3d'
%
% OUTPUT
%    s       : a string in TeX-Notation representing the number in
%              sientific notation
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-03-29
% see LICENSE.txt

    expn = floor( log( abs(x) )/log(10) ); % prediction
    expn = expn + ( x == 10^(expn+1) ); % correction
    mant = x / 10^expn;
    if (expn~=0)
        s = sprintf([fmtmant,'\\times 10^{',fmtexpn,'}'], mant, expn);
    else
        s = sprintf(fmtmant, mant); 
    end
end

