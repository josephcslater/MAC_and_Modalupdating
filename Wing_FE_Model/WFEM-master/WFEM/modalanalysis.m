function [fs,fms]=modalanalysis
% MODALANALYSIS solves for the modes of the structure and plots
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
lslave=length(slavedofs);

if isempty(fs)|restart==0,% This is old stuff from previous restart
                          % code. If we haven't any natural
                          % frequencies already calculated, or we
                          % don't want a restart.
  if lslave>0% Here we need to reduce out constrained DOFs before
             % static analysis. Slave DOFs are defined when we apply
             % constraints or rigid elements.
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
    [Mr,Kr,Tr,master,slave]=guyan(M,K,master);
  else
    [Mr,Kr,Tr,master,slave]=guyan(M,K);
  end
  % From $$$$$$ to $$$$$$ is the actual eigensolution.
  OPTS.issym=1;
  OPTS.isreal=1;
  %To increase accuracy: (Agnes)
  %OPTS.tol=eps/10;
  [Tr,Mr,Kr]=reducedofs;
  [minvals,minvallocs]=sort(diag(Kr)./diag(Mr));
  
  shift=minvals(min(7,length(minvals)));
  [fms,f]=eigs((Kr+Kr')/2+shift*(Mr+Mr')/2,(Mr+Mr')/2,min([ size(Kr,1) ...
		    max([floor(sqrt(size(Kr,1))) 100])]),0,OPTS);
  fs=sqrt(diag(f)-shift)/2/pi;
  % $$$$$$
  fms=Tr*fms;
  [fs,fsindex]=sort(fs);
  fms=fms(:,fsindex);
else
  disp('Using previous eigensolution.')
end


