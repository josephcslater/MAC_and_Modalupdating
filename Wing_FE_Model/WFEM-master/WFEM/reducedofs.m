function [Tr,Mr,Kr,Fnlnr,Fr]=reducedofs
% [Tr,Mr,Kr,Fnlnr,Fr]=REDUCEDOFS reduces the number of dofs
% Also will perform all force projections on reduced DOFs if
% wanted. 
  
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

i=1:size(K,1);
lslave=length(slavedofs);
lzeros=sum(diag(K)==0);

if lslave>0|lzeros>0
  i(slavedofs)=zeros(length(slavedofs),1);%pop zeros into the k
                                          %indices corresponding
                                          %to slave dofs
  i=sort(i); % sort them. This puts the slave coordinates all up
             % top. 
  master=i(lslave+1:size(nodes,1)*6);%The master coordinates are
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
  [Mr,Kr,Tr,master,slave,zerodofs]=guyan(M,K,master);
  
    
else
  Mr=M;
  Kr=K;
  Tr=eye(size(M));
  master=(1:length(M))';
end
