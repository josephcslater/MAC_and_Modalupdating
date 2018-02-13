function out=panel1(mode,b,c,d,e)

% PANEL1 denotes that panel has some particular geometrical
% interest. This is useful for evaluating surface quality for
% radar-like applications.
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
global nout
global curlineno

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
  out=4;
end
if strcmp(mode,'generate')
%  This just adds this element to the element data structure. 
  elnum=c;
  
  element(elnum).nodes=[b(1) b(2) b(3) b(4)];
  
  element(elnum).properties=b(5);
  element(elnum).lineno=curlineno;
  %How should I draw the darn thing. Likely a surface. Right. Using
  %the draw call, or not?
  panelcolor=[1 0 1];
  %Don't like this color? Use colorui to pick another one. Another
  %option is that if we can't see the elements separately, we can
  %chunk up x*y*z, divide by x*y*x of element, see if we get
  %integer powers or not to define colors that vary by panel. 
  surfs=[surfs;element(elnum).nodes panelcolor];
end
if strcmp(mode,'make')
  elnum=b;
  elnodes=element(elnum).nodes;
  elmass=elprops(element(elnum).properties).a;
  elnodalmass=elmass/4;
  %elmassmatrix=speye(12)*elnodalmass;
  elmassmatrix=diag([1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0])*elnodalmass;
  indices=[elnodes(1)*6+(-5:0) elnodes(2)*6+(-5:0) elnodes(3)*6+(-5:0) ...
	   elnodes(4)*6+(-5:0) ];
  M(indices,indices)=M(indices,indices)+elmassmatrix;
  % There is nothing to make with a surface node element. Its
  % existence is everything it does.
  nout=size(Cd,1);
  Cd(nout+(1:6),elnodes(1)*6+(-5:0))=eye(6);
  Cd(nout+6+(1:6),elnodes(2)*6+(-5:0))=eye(6);
  Cd(nout+12+(1:6),elnodes(3)*6+(-5:0))=eye(6);
  Cd(nout+18+(1:6),elnodes(4)*6+(-5:0))=eye(6);

  
  %   Xd(scd+(1:24),1)=[...
  %       nodes(elnodes(1),1);...
  %       nodes(elnodes(1),2);...
  %       nodes(elnodes(1),3);... 
  %       nodes(elnodes(2),1);...
  %       nodes(elnodes(2),2);...
  %       nodes(elnodes(2),3);... 
  %       nodes(elnodes(3),1);...
  %       nodes(elnodes(3),2);...
  %       nodes(elnodes(3),3);... 
  %       nodes(elnodes(4),1);...
  %       nodes(elnodes(4),2);...
  %       nodes(elnodes(4),3)];
  
  Xd(nout+(1:24),1)=[nodes(elnodes(1),:)';0;0;0; nodes(elnodes(2),:)';0;0;0; nodes(elnodes(3),:)';0;0;0; nodes(elnodes(4),:)';0;0;0 ];
  nout=nout+24;
  
  
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
