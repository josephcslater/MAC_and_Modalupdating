function [Ass, Bss, Css, Dss, Z0, Tr]=genssm(modenums)
% [Ass, Bss, Css, Dss, Z0,Tr]=genssm(MODENUMS) generates a state space
% model. Tr post-multiplied by displacements or velocities yields
% the complete nodel displacement vector or complete nodal velocity
% vector. In order to extract portions of the motion
% (e.g. non-rigid body portions) one should simulate with an
% identity matric for C (note that this code uses the variable 'C'
% for the damping matrix, Css is the state-space output matrix)),
% followed by removal of the uninteresting parts, then multiply by
% the generated C matrix. 
% MODENUMS are the modes to keep when doing modal reduction (see
% modalreduction). The default is 7:30 (truncating the first 6
% modes, an all modes higher than 30). 
  
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
global Tr  Mr Kr
global C
global X0 V0
global Cd Cv Ca Bso Xd
global alpha beta
global nout nin


i=1:size(K,1);
lslave=length(slavedofs);

if lslave>0% Here we need to reduce out constrained DOFs before
	   % generating the state space model. Slave DOFs are defined when we apply
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
disp('entering modalreduction')
if nargin==0
  modenums=7:30;
end
if modenums(1)==0
    modenums=modenums(2:length(modenums));
end


[Mred,Kred,T]=modalreduction(Mr,Kr,modenums);
disp('Model Modallly Reduced')

%
Tr=Tr*T;
Mr=Mred;
Kr=Kred;

sm=size(Mr,1);
%%%%%%%%%%%%%%%%%
% Need to make a damping matrix here using proportional
% damping. Uncomment the next three lines. Once I define a
% proportional damping element, only the last line will be used.
%alpha=.1;
%beta=.1;
%C=alpha*M+beta*K;


if exist('C')==1
  Cr=Tr'*C*Tr;
else
  Cr=spalloc(sm,sm,0);
end
%more on
%[eig(full(Mr)) eig(full(Kr)) eig(full(Cr))]
%save temp2
%pause
Ass=[spalloc(sm,sm,0) speye(sm);...
   -Mr\Kr   -Mr\Cr];

% Just in case nothing was ever assigned to an input or output, we
% have to make sure the matrices are the right size.
if size(Ca,1)<nout
  Ca(nout,1)=0;
end
if size(Cd,1)<nout
  Cd(nout,1)=0;
end
if size(Cv,1)<nout
  Cv(nout,1)=0;    
end

if size(Ca,2)<nin
%  if nin>0
%    Ca(1,nin)=0;
%  else
    Ca=sparse(nin,nout);
%  end
end
Bso(size(Tr,1),1)=0;
%size(nodes,1)*6
%disp('a')


%size(Tr')
%size(Bso)
Bsor=Tr'*Bso;
Bss=[spalloc(size(Bsor,1),size(Bsor,2),0);
     Mr\Bsor];
Ca(1,size(Tr,1))=0;
Cv(1,size(Tr,1))=0;
Cd(1,size(Tr,1))=0;disp('sizes')
Car=Ca*Tr;size(Car);
Cvr=Cv*Tr;size(Cvr);
Cdr=Cd*Tr;size(Cdr);
size(Mr);size(Kr);

Css=[Cdr-Car/Mr*Kr Cvr-Car/Mr*Cr];
disp('Min eig of Mr')
min(eig(full(Mr)))
Dss=Car/Mr*Bsor;
X0(size(Tr,1))=0;
V0(size(Tr,1))=0;
%Z0=[Tr'*X0;Tr'*V0];

Z0=[];
