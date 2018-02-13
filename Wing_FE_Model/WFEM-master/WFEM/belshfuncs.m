% Derive Shape functions for quadratic beam element
% For all points, the polynomials N and N' are evaluated.
% w=a+bx+cx^2+dx^3+ex^4+fx^5
b=[1 -1  1 -1  1 -1; % w(-1) x^n
   1  1  1  1  1  1; % w(0) x^n
   1  0  0  0  0  0; % w(1) x^n
   0  1 -2  3 -4  5; % w'(-1) x^n
   0  1  2  3  4  5; % w'(0) x^n
   0  1  0  0  0  0];% w'(1) x^n 
b=fliplr(b);
disp('b flipped to decreasing power order')
% For N1
a=[1 0 0 0 0 0]';
c1=(b\a)'
c1d=polyderiv(c1)
c1dd=polyderiv(polyderiv(c1))
% For N2
a=[0 0 0 1 0 0]';
c2=(b\a)'
c2d=polyderiv(c2)
c2dd=polyderiv(polyderiv(c2))
% For N2
a=[0 1 0 0 0 0]';
c3=(b\a)'
c3d=polyderiv(c3)
c3dd=polyderiv(polyderiv(c3))
% For N4
a=[0 0 0 0 1 0]';
c4=(b\a)'
c4d=polyderiv(c4)
c4dd=polyderiv(polyderiv(c4))
% For N6
a=[0 0 1 0 0 0]';
c5=(b\a)'
c5d=polyderiv(c5)
c5dd=polyderiv(polyderiv(c5))
% For N5
a=[0 0 0 0 0 1]';
c6=(b\a)'
c6d=polyderiv(c6)
c6dd=polyderiv(polyderiv(c6))
i=-1:.01:1;
plot(i,polyval(c1,i),i,polyval(c2,i),i,polyval(c3,i),...
    i,polyval(c4,i),i,polyval(c5,i),i,polyval(c6,i))

%plot(i,polyval(c1,i),"-;sf1;",i,polyval(c2,i),"-;sf2;",i,polyval(c3,i),"-;sf3;",i,polyval(c4,i),"-;sf4;",i,polyval(c5,i),"-;sf5;",i,polyval(c6,i),"-;sf6;")
