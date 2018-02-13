function J = jacobian(X,Y,Z,pt)
% Evaluate the Jacobian at point pt (x,y,z)
% X, Y, and Z are the coordinates of 8 nodes in the following form
% X = [x1, x2, x3, ..., x8]
% Y = [y1, y2, y3, ..., y8]
% Z = [z1, z2, z3, ..., z8]

% Shape Function and corresponding derivative
Results = brick8_ShapeFun();

bn1 = Results.c1;
bn1d = Results.c1d;

bn2 = Results.c2;
bn2d = Results.c2d; 

bn3 = Results.c3;
bn3d = Results.c3d;

bn4 = Results.c4;
bn4d = Results.c4d;  

bn5 = Results.c5;
bn5d = Results.c5d;

bn6 = Results.c6;
bn6d = Results.c6d;

bn7 = Results.c7;
bn7d = Results.c7d;

bn8 = Results.c8;
bn8d = Results.c8d; 

% Jacobian components in the 3D plynomial form 
% w=a1+a2x+a3y+a4z+a5xy+a6xz+a7yz+a8xyz
Jp11 = bn1d(1,:)*X(1) + bn2d(1,:)*X(2) + bn3d(1,:)*X(3) + bn4d(1,:)*X(4) + ...
      bn5d(1,:)*X(5) + bn6d(1,:)*X(6) + bn7d(1,:)*X(7) + bn8d(1,:)*X(8);
Jp12 = bn1d(1,:)*Y(1) + bn2d(1,:)*Y(2) + bn3d(1,:)*Y(3) + bn4d(1,:)*Y(4) + ...
      bn5d(1,:)*Y(5) + bn6d(1,:)*Y(6) + bn7d(1,:)*Y(7) + bn8d(1,:)*Y(8);
Jp13 = bn1d(1,:)*Z(1) + bn2d(1,:)*Z(2) + bn3d(1,:)*Z(3) + bn4d(1,:)*Z(4) + ...
      bn5d(1,:)*Z(5) + bn6d(1,:)*Z(6) + bn7d(1,:)*Z(7) + bn8d(1,:)*Z(8); 
  
Jp21 = bn1d(2,:)*X(1) + bn2d(2,:)*X(2) + bn3d(2,:)*X(3) + bn4d(2,:)*X(4) + ...
      bn5d(2,:)*X(5) + bn6d(2,:)*X(6) + bn7d(2,:)*X(7) + bn8d(2,:)*X(8);
Jp22 = bn1d(2,:)*Y(1) + bn2d(2,:)*Y(2) + bn3d(2,:)*Y(3) + bn4d(2,:)*Y(4) + ...
      bn5d(2,:)*Y(5) + bn6d(2,:)*Y(6) + bn7d(2,:)*Y(7) + bn8d(2,:)*Y(8);
Jp23 = bn1d(2,:)*Z(1) + bn2d(2,:)*Z(2) + bn3d(2,:)*Z(3) + bn4d(2,:)*Z(4) + ...
      bn5d(2,:)*Z(5) + bn6d(2,:)*Z(6) + bn7d(2,:)*Z(7) + bn8d(2,:)*Z(8);
  
Jp31 = bn1d(3,:)*X(1) + bn2d(3,:)*X(2) + bn3d(3,:)*X(3) + bn4d(3,:)*X(4) + ...
      bn5d(3,:)*X(5) + bn6d(3,:)*X(6) + bn7d(3,:)*X(7) + bn8d(3,:)*X(8);  
Jp32 = bn1d(3,:)*Y(1) + bn2d(3,:)*Y(2) + bn3d(3,:)*Y(3) + bn4d(3,:)*Y(4) + ...
      bn5d(3,:)*Y(5) + bn6d(3,:)*Y(6) + bn7d(3,:)*Y(7) + bn8d(3,:)*Y(8);   
Jp33 = bn1d(3,:)*Z(1) + bn2d(3,:)*Z(2) + bn3d(3,:)*Z(3) + bn4d(3,:)*Z(4) + ...
      bn5d(3,:)*Z(5) + bn6d(3,:)*Z(6) + bn7d(3,:)*Z(7) + bn8d(3,:)*Z(8); 

% Evaluate the Jacobian at point pt
J11 = poly3dval(Jp11,pt);
J12 = poly3dval(Jp12,pt);
J13 = poly3dval(Jp13,pt);

J21 = poly3dval(Jp21,pt);
J22 = poly3dval(Jp22,pt);
J23 = poly3dval(Jp23,pt);

J31 = poly3dval(Jp31,pt);
J32 = poly3dval(Jp32,pt);
J33 = poly3dval(Jp33,pt);
  
J = [J11, J12, J13;
     J21, J22, J23;
     J31, J32, J33];
  
  
  
  


