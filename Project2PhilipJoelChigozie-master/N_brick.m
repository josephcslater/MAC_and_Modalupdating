function [ N ] = N_brick( Epsilon, Eta, Zeta )
% N_brick uses inputs Epsilon, Eta and Zeta to generate N

% Epsilon is Xi in the notes. EpsilonI, EtaI, and ZetaI are the nodal
% values in order
EpsilonI=[-1 -1 -1 -1 1 1 1 1];
EtaI=[-1 -1 1 1 -1 -1 1 1];
ZetaI=[1 -1 -1 1 1 -1 -1 1];
N=[]; % creates N matrix
for i=1:8 % iterates over each node
    Ni=(1/8)*(1+Epsilon*EpsilonI(i))*(1+Eta*EtaI(i))*(1+Zeta*ZetaI(i)); % Calculates each nodal shape function
    Ne=[Ni Ni Ni]; % Creates a nodal N matrix
    Ne=diag(Ne); % stores the shape function on the diagonal
    N=[N Ne]; % appends the element shape function matrix with the nodal matrix
end

end

