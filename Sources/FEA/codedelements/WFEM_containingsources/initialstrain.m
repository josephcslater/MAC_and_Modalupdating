function X=initialstrain
% INITIALSTRAIN solves the static response problem and outputs the
% results to a file. Right now, those results are nothing. We'll
% add to that later. Next, the results are plotted to the
% screen. Why do we need a whole code for this? Well, most of this
% routine will be defining how to output the data to a file.
% This routine will also be a precurser to the stress/strain
% operations for the static solution. 
  
global K
global Ks
global M
global nodes
global beamprops
global matprops
global points
global elprops
global element
global Fepsn % Forces representing initial strain
global F     % Prescribed forces
global Fnln  % Forces due to changing coordinate systems. 
global lines
global surfs
global bcs %Boundary Conditions
size(F);
size(Fepsn);%size(K)
%[v,d]=eigs(K,size(K,1));%size(d)

%The way it should be done with rigid body modes
% [v,d]=eig(full(K));
% diagd=diag(d);
% id=diag(1./diagd.*(abs(diagd)>1000*eps));
% X=v*id*v'*(Fepsn);


%Below works if K not singular
X=K\Fepsn;

%Truncate Lagrange Multiplier 'deflections'
X=X(1:size(nodes,1)*6,1);

%Dump any imaginary part
X=real(X);
%plotgeom('deformed',nodes,lines,X)

%figure(gcf)