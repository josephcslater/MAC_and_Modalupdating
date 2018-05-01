function remove(removeddofs)

% REMOVE removes the DOFS requested from the FEM. It does this the
% elementary way, by removing rows and columns. It notes the
% removal by mapping all removed DOFS to a model dof higher than
% any real dof. Assignment of forces after initial model creation
% must used the dofs vector to look up the new location of
% DOFS. That is, the ith element of the variable DOF is the
% internal DOF number after model reduction. 
% This is really a bad idea, and I would recomment using the
% CONSTRAIN code instead. Only problem being it dosn't exists at
% the time of this writing.
  
% Copyright Joseph C. Slater, 2002
  
  
global dofs
global ndofs
global K
global M
dofs
dofs(:,2)=zeros(ndofs,1);
dofs(removeddofs,2)=ones(length(removeddofs),1);
[s,i]=sort(dofs(:,2));
dofs=dofs(i,1);
retaineddofs=dofs(1:(ndofs-length(removeddofs)));
K=K(retaineddofs,retaineddofs);
M=M(retaineddofs,retaineddofs);