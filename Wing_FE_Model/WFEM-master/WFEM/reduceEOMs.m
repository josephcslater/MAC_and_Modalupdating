function [Mr,Kr,Tr]=reduceEOMs
% [Mr,Kr,Tr]=reduceEOMs removes constraint coordinates AND removes zero stiffness DOFs
% from truss elements, etc. 
  
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
global slavedofs

size(F);
size(Fepsn);

i=1:size(K,1);
lslave=length(slavedofs);
  
if lslave>0% Here we need to reduce out constrained DOFs before
    % eigen analysis. Slave DOFs are defined when we apply constraints or
    % rigid elements.
    i(slavedofs)=zeros(length(slavedofs),1);%pop zeros into the k
    %indices corresponding to slave dofs
    i=sort(i); % sort them. This puts the slave coordinates all up
    % top.
    master=i(lslave+1:size(nodes,1)*6);%The master coordinates are
    %the ones that we didn't reduce out via some constraints. The master
    %coordinates are all of the coordinates after the zero values in i, but
    %not those beyond coordinates corresponding to nodal DOFS. The upper
    %limit of size(nodes,1)*6 fixes this at the last real DOF.
    [Mr,Kr,Tr,master,slave]=guyan(M,K,master);
else
    [Mr,Kr,Tr,master,slave]=guyan(M,K);
end
%[Mr,Kr,Tr,master,slave]=guyan(M,K,master);
% Fr=Tr'*(F+Fepsn);
% Xr=Kr\Fr;
% X=Tr*Xr;
