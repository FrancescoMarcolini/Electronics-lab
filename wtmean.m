function [xm,dxm,chi2] = wtmean(x,dx)
% [xm,dxm] = wtmean(x,dx)
%
% Computes the mean of x weighted by the uncertanties dx.
%
% If called in the form [xm,dxm,chi2] = wtmean(x,dx) it computes also che
% chi square divided by the number of degrees of fredom.

w = 1./dx.^2;
sw = sum(w);
sxw = sum(x.*w);
xm = sxw/sw;
dxm = sqrt(1/sw);

if nargout == 3
    N = length(x);
    chi2 = sum(((x - xm).^2).*w)/(N-1);
end