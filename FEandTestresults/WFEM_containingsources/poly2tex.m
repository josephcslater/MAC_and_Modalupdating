function texstr=poly2tex(poly,var,iseqn)
% POLY2TEX(POLY,VAR,ISEQN) returns the latex equation for a polynomial
% POLY is the polynomial defined in matlab/octave form. 
% VAR is the variable the polynomial is to be written in.
% ISEQN (Optional) returns the polynomial equal to zero (an
% equation) instead of just the polynomial.
  
% Copyright Joseph C. Slater, 2002
if nargin==2
  iseqn=0;
end
n=length(poly)-1;
texstr='';
if n>1
  for i=n:-1:2
    texstr=[texstr num2str(poly(n-i+1)) ' ' var '^{' num2str(i) ...
                        '} + '];
  end
end

texstr=[texstr  num2str(poly(n)) ' ' var  ' + '];
texstr=[texstr  num2str(poly(n+1)) ];
if iseqn==1
  texstr=[texstr ' = 0'];
end

