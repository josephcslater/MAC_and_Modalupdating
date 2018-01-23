function [ Emat ] = E_matrix( E,G )
% E_matrix solves for Emat given E and G
nu=(E/(2*G))-1; % Poisson's Ratio
val1=1-nu; % Calculates values used in the matrix
val2=1+nu;
val3=1-2*nu;
Emat=zeros(6,6); % Creates a matrix of all zeros
Emat(1:3,1:3)=nu; % Fills the top left 3x3 block with nu
Emat(1,1)=val1; % Overwrites the diagonal
Emat(2,2)=val1;
Emat(3,3)=val1;
Emat(4,4)=val3/2; % Writes the bottom right diagonal
Emat(5,5)=val3/2; 
Emat(6,6)=val3/2;
Emat=Emat*E/(val2*val3); % Multiplies by the outside constants

end

