function y = poly3dval(p,X)
% Evaluate the value of the following polynomial
% y = p1 + p2*x + p3*y + p4*z + p5*xy + p6*xz + p7*yz + p8*xyz
% with x, y and z as
% x = X(1)
% y = X(2)
% z = X(3)

x = X(1);
y = X(2);
z = X(3);

y = p(1) + p(2)*x + p(3)*y + p(4)*z + p(5)*x*y + p(6)*x*z + p(7)*y*z + p(8)*x*y*z;



