function [modes]=buildrbmodes
% [modes]=buildrbmodes builds the rigid body modes 'cleanly', so
% that they are pure rotations about axes and pure
% translations. Further, they are normalized. They are NOT in
% reduced coordinates. 
  
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


%[totalmass,inertiatensor, cg]=findinertia;
snodes=size(nodes,1);
emptymode=zeros(snodes,1);
%Mode 1 is diaplacement in the x direction
mode1=emptymode;
mode1(1:6:snodes*6)=ones(snodes,1)/sqrt(snodes);
disp('Mode 1 built')
%Mode 2 is diaplacement in the y direction
mode2=emptymode;
mode2(2:6:snodes*6)=ones(snodes,1)/sqrt(snodes);
disp('Mode 2 built')
%Mode 3 is diaplacement in the z direction
mode3=emptymode;
mode3(3:6:snodes*6)=ones(snodes,1)/sqrt(snodes);
disp('Mode 3 built')
%Mode 4 is rotation about the x axis
mode4=emptymode;

  %rotational coordinates
  mode4(4:6:snodes*6)=ones(snodes,1)/sqrt(snodes); 
  
  %y displacements;
  mode4(2:6,snodes*6)=-nodes(:,3);
  
  %z displacements;
  mode4(3:6,snodes*6)=nodes(:,2);
  
  mode4=mode4/norm(mode4);
  
disp('Mode 4 built')
%Mode 5 is rotation about the y axis
mode5=emptymode;

  %rotational coordinates
  mode5(5:6:snodes*6)=ones(snodes,1)/sqrt(snodes); 
  
  %x displacements;
  mode5(1:6,snodes*6)=nodes(:,3);
  
  %z displacements;
  mode5(3:6,snodes*6)=-nodes(:,1);
  
  mode5=mode5/norm(mode5);
  
  disp('Mode 5 built')
%Mode 6 is rotation about the z axis
mode6=emptymode;

  %rotational coordinates
  mode6(6:6:snodes*6)=ones(snodes,1)/sqrt(snodes); 
  
  %y displacements;
  mode6(2:6,snodes*6)=nodes(:,1);
  
  %x displacements;
  mode6(1:6,snodes*6)=-nodes(:,2);
  
  mode6=mode6/norm(mode6);
  disp('Mode 6 built')
