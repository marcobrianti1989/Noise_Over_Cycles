%
%
% c = cfilter(x,P0,P1) is the bandpass filter extracting from the 
% series x the component with periodicity [P0,P1] and stores the result
% in the column vector y. x must be a column vector, while P0 and
% P1 must be real (scalar) numbers such that P0<P1. If n is the number
% of entries in x, then y has n-1 entries (the first observation
% is lost). 
% [c, T] = cfilter(x,P0,P1) produces also T = x - c.
%
%
function [c,T] = cfilter(x,P0,P1)
%
%
% compute the first difference of x and call it dc.
%
dx = diff(x)-mean(diff(x));
%
% compute the limits of the frequency band, i.e.
% the scalars lambda0 and lambda1 
%
lambda0 = 2*pi/P1;
lambda1 = 2*pi/P0;
%
%
% compute the vector h, whose elements are the cofficients
% of the C-filter. The theoretical filter is truncated at lags n
% and -n.
%
n = length(x);
w = (sin((1:n)*lambda1)'-sin((1:n)*lambda0)')./(1:n)';
h = [flipud(w/pi);(lambda1-lambda0)/pi;w/pi];
%
%
% apply the filter to the series dx to get the 
% first difference of the cycle dc. dc has the same 
% length as dx and is obtained by truncating the last n 
% entries from both sides of the convolution of h and dx. 
%
yy = conv(h,dx);
dc = yy(n+1:length(yy)-n);
%
%
% compute c as the centered cumulated sum of dc and T 
% as the difference between the last n-1 observations of x and  c.
%
c = cumsum(dc);
c = c-mean(c); 
T = x(2:length(x))-c;