function [totalmass,inertiatensor, cg]=findinertia
% [totalmass,INERTIATENSOR,CG]=FINDINERTIA solves for the inertia tensor
% relative to the center of gravity, as well as the center of gravity
% them to the screen with a pretty lame menu. More later.
  
global restart
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
global bcs   % Boundary Conditions
global slavedofs
global fs
global fms 
global filename
global Tr 
global Mr
global Kr



i=1:size(K,1);
%Reduce the mass and stiffness matrices via Guyan. For a system
%with 6 rigid body modes, the stiffness matrix shoul be all zeros. 
lslave=length(slavedofs);
slavedofs;
if lslave>0
  % Here we need to reduce out constrained DOFs before
  % static analysis. Slave DOFs are defined when we apply
  % constraints or rigid elements.
  
  i(slavedofs)=zeros(length(slavedofs),1);
  %pop zeros into the k
  %indices corresponding
  %to slave dofs
  i=sort(i); 
  % sort them. This puts the slave coordinates all up
  % top. 
  master=i(lslave+1:size(nodes,1)*6);
  %The master coordinates are
  %the ones that we didn't
  %reduce out via some
  %constraints. The master
  %coordinates are all of the
  %coordinates after the zero
  %values in i, but not those
  %beyond coordinates
  %corresponding to nodal
  %DOFS. The upper limit of
  %size(nodes,1)*6 fixes this
  %at the last real DOF.
  
  
  %master=7:12;
  [Mr,Kr,Tr,master,slave]=guyan(M,K,master);
  
else
  [Mr,Kr,Tr,master,slave]=guyan(M,K);
end

full(Mr);

%Find the total mass
mnodevals=sum(M(1:6:size(nodes,1)*6,1:6:size(nodes,1)*6))';
totalmass=sum(mnodevals);
mnodevals2=sum(M(2:6:size(nodes,1)*6,2:6:size(nodes,1)*6))';
totalmass=sum(mnodevals);
mnodevals3=sum(M(3:6:size(nodes,1)*6,3:6:size(nodes,1)*6))';
totalmass=sum(mnodevals);

% Now to find the CG. 


xbar=(nodes(:,1)'*mnodevals)/totalmass;
ybar=(nodes(:,2)'*mnodevals)/totalmass;
zbar=(nodes(:,3)'*mnodevals)/totalmass;


M=Mr(1:3,1:3);M2=Mr(4:6,1:3);I=Mr(4:6,4:6);

inertiatensor=full(I+M2/M*M2);


%I=Mr(4:6,4:6)-diag([xbar,ybar,zbar])*totalmass*ones(3,3)*diag([xbar,ybar,zbar]);

%inertiatensor=I;
cg=[xbar,ybar,zbar];


