function dp = polyderiv(p)
% Take one derivative of a polynomial
% Copyright Joseph C. Slater
% Oct. 8, 2016

l = length(p);
pows = (l-1):-1:0;
dp = pows.*p;
dp = dp(1:(length(dp)-1));