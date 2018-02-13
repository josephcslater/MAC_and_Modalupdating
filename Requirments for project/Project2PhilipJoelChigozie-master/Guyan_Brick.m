function [ Kr,Mr ] = Guyan_Brick( Ke,Me )
% Guyan_Brick reduces the stiffness and mass matrix size to 24x24 using
% guyan reduction

% check to see if Ke needs to be reduced
if size(Ke,1)>24 || size(Ke,2)>24
    K11=Ke(1:24,1:24);
    K21=Ke(25:end,1:24);
    K22=Ke(25:end,25:end);
    T=[eye(size(K11,1));-K22\K21];
    Kr=T'*Ke*T;
else
    Kr=Ke;
end

% check to see if Me needs to be reduced
if size(Me,1)>24 || size(Me,2)>24
    M11=Me(1:24,1:24);
    M21=Ke(25:end,1:24);
    M22=Me(25:end,25:end);
    T=[eye(size(M11,1));-M22\M21];
    Mr=T'*Me*T;
else
    Mr=Me;
end
end

