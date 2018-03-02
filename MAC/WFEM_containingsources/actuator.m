function out=actuator(mode,b,c,d,e)

% Actuator places an input at the defined location, in the defined
% direction.
% Format:
% actuator elements
% node# direction
%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables (global):
% -------------------
% K       :    Global stiffness matrix
% Ks      :    Global stiffness buckling matrix
% M       :    Global mass matrix
% nodes   :    [x y z] nodal locations

global K
global Ks
global M
global nodes % Node locations
global elprops
global element
global points
global surfs
global Fepsn % Initial strain "forces". 
global lines
global Cd Cv Ca Xd
global nin
global curlineno
global Bso %Second order form input matrix

%
% Variables (local):
% ------------------
% Cripe, the code is only 50 lines long. I think we can spare this list.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Joseph C. Slater, 10/13/2003.
% joseph.slater@wright.edu
out=0;
if strcmp(mode,'numofnodes')
  out=1;
end
if strcmp(mode,'generate')
%  This just adds this element to the element data structure. 
  elnum=c;
  element(elnum).nodes=[b(1)];% b(2) b(3) b(4)];
  
  element(elnum).properties=b(2);
%  element(elnum).lineno=curlineno;
  %How should I draw the darn thing. 
end
if strcmp(mode,'make')
  elnum=b;
  elnodes=element(elnum).nodes;
  %element(elnum)
  eldirection=element(elnum).properties;

  %indices=[elnodes(1)*6+(-5:0) elnodes(2)*6+(-5:0) elnodes(3)*6+(-5:0) ...
	   %elnodes(4)*6+(-5:0) ];
  
  %sbs=size(Bso,2);
  nin=nin+1;
  %elnodes
  %eldirection
  %elnodes
  %eldirection
  
  Bso(((elnodes-1)*6)+eldirection,nin)=1;
  
  
  %size(Bso)
  %pause
% $$$   Cd(scd+(1:6),elnodes(1)*6+(-5:0))=eye(6);
% $$$   Cd(scd+6+(1:6),elnodes(2)*6+(-5:0))=eye(6);
% $$$   Cd(scd+12+(1:6),elnodes(3)*6+(-5:0))=eye(6);
% $$$   Cd(scd+18+(1:6),elnodes(4)*6+(-5:0))=eye(6);

  
  
  %Xd(scd+(1:24),1)=[nodes(elnodes(1),:)';0;0;0; nodes(elnodes(2),:)';0;0;0; nodes(elnodes(3),:)';0;0;0; nodes(elnodes(4),:)';0;0;0 ];
      
  %  Cv(scd+(1:6),elnodes(1)*6+(-5:0))=zeros(6);
  %  Cv(scd+6+(1:6),elnodes(2)*6+(-5:0))=zeros(6);
  %  Cv(scd+12+(1:6),elnodes(3)*6+(-5:0))=zeros(6);
  %  Cv(scd+18+(1:6),elnodes(4)*6+(-5:0))=zeros(6);
  %  if size(Ca,1)==0
  %    Ca=sparse(24,size(Ca,2));
  %  end    
  %  Ca(size(Ca,2)
elseif strcmp(mode,'draw')
  elnodes=element(elnum).nodes;
  %What benefit to changing colors? Use generic?
elseif strcmp(mode,'buckle')
end
