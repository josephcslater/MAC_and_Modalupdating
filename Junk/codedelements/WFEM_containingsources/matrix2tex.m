function matrix2tex(a,format)
% MATRIX2TEX(MATRIX,PLACES) returns the latex formatted matrix
% MATRIX is the polynomial defined in matlab/octave form.
% Format is a string representing the fprint format
% %d Decimal notation (e.g. integer)
% %e Exponential notation (% in 3.1415e+00)
% %E Exponential notation (using an uppercase E as in 3.1415E+00)
% %f Fixed-point notation
% %g The more compact of %e or %f. Insignificant zeros do not print, sorry.
% %G Same as %g, but using an uppercase E
% etc. See C reference for more information.
% Numbers (two number to be exact, with a decimal between them) can be injected between
% the '%' and the letter, as per a C reference of your choice.
% Default format is '%5.3f';
%
% Example:
% a=rand(4,2);
% matrix2tex(a)
% matrix2tex(a,'%8.3g')
% a=a*1000;
% matrix2tex(a,'%8.3e')

% Copyright Joseph C. Slater, 2005
% Updated 7/10/09 to handle complex numbers

texstr='';
m=size(a,1);
n=size(a,2);
if nargin==1
	format='%5.3f';
end

eol='\\\\';
formataug=[format '&'];


for i=1:m
	for j=1:n
		if real(a(i,j))~=0|imag(a(i,j))==0
			texstr=[texstr  sprintf(format,real(a(i,j))) ];
		end
		if sign(imag(a(i,j)))==1&real(a(i,j))~=0&imag(a(i,j))~=0
			texstr=[texstr '+'];
		end
		if imag(a(i,j))~=0
			texstr=[texstr  sprintf(format,imag(a(i,j))) 'i'  ];
		end
		if j==n&i==m
			texstr=[texstr  ' \n'];
		elseif j==n
			texstr=[texstr  eol ' \n'];
		else
			texstr=[texstr '&'];
		end
	end
end




%Octave/matlab differentiator
if exist('ver')==6
	sprintf(texstr)
else
	printf(texstr)
end




