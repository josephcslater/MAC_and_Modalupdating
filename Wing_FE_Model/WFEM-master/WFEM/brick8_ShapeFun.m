function Results = brick8_ShapeFun()
% Derive Shape functions for brick8 element
% For all points, the polynomials N are evaluated.
% w=a1+a2x+a3y+a4z+a5xy+a6xz+a7yz+a8xyz
b=[ 1,-1,-1, 1, 1,-1,-1, 1; % w(-1,-1, 1) @ point 1
    1,-1,-1,-1, 1, 1, 1,-1; % w(-1,-1,-1) @ point 2 
    1,-1, 1,-1,-1, 1,-1, 1; % w(-1, 1,-1) @ point 3
    1,-1, 1, 1,-1,-1, 1,-1; % w(-1, 1, 1) @ point 4
    1, 1,-1, 1,-1, 1,-1,-1; % w( 1,-1, 1) @ point 5
    1, 1,-1,-1,-1,-1, 1, 1; % w( 1,-1,-1) @ point 6
    1, 1, 1,-1, 1,-1,-1,-1; % w( 1, 1,-1) @ point 7
    1, 1, 1, 1, 1, 1, 1, 1]; % w( 1, 1, 1) @ point 8
                     

% b=fliplr(b);
% disp('b flipped to decreasing power order')
% For N1
a=[1 0 0 0 0 0 0 0]';
c1=(b\a)';
c1d=brick8_polyderiv(c1);

% For N2
a=[0 1 0 0 0 0 0 0]';
c2=(b\a)';
c2d=brick8_polyderiv(c2);
% For N3
a=[0 0 1 0 0 0 0 0]';
c3=(b\a)';
c3d=brick8_polyderiv(c3);
% For N4
a=[0 0 0 1 0 0 0 0]';
c4=(b\a)';
c4d=brick8_polyderiv(c4);
% For N5
a=[0 0 0 0 1 0 0 0]';
c5=(b\a)';
c5d=brick8_polyderiv(c5);
% For N6
a=[0 0 0 0 0 1 0 0]';
c6=(b\a)';
c6d=brick8_polyderiv(c6);
% For N7
a=[0 0 0 0 0 0 1 0]';
c7=(b\a)';
c7d=brick8_polyderiv(c7);
% For N8
a=[0 0 0 0 0 0 0 1]';
c8=(b\a)';
c8d=brick8_polyderiv(c8);

Results = struct('c1',c1,'c1d',c1d,...
                 'c2',c2,'c2d',c2d,...
                 'c3',c3,'c3d',c3d,...
                 'c4',c4,'c4d',c4d,...
                 'c5',c5,'c5d',c5d,...
                 'c6',c6,'c6d',c6d,...
                 'c7',c7,'c7d',c7d,...
                 'c8',c8,'c8d',c8d);


function dp = brick8_polyderiv(a)
% Take derivative of the following polynomial with respect to x, y and z
% w=a1+a2x+a3y+a4z+a5xy+a6xz+a7yz+a8xyz
dp_x = [a(2), 0, a(5), a(6), 0, 0, a(8), 0];
dp_y = [a(3), a(5), 0, a(7), 0, a(8), 0, 0];
dp_z = [a(4), a(6), a(7), 0, a(8), 0, 0, 0];
dp = [dp_x; dp_y; dp_z];



