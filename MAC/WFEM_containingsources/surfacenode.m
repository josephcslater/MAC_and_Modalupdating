function out=surfacenode(mode,b,c,d,e)
  
% SURFACENODE denotes that this node is on some surface. thi is
% useful for evaluating surface quality for radar-like aplicatoins.
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
global Fepsn % Initial strain "forces". 
global lines

%
% Variables (local):
% ------------------
% Cripe, the code is only 50 lines long. I think we can spare this list.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Joseph C. Slater, 10/15/2002.
% joseph.slater@wright.edu
out=0;
if strcmp(mode,'numofnodes')
  out=1;
end
if strcmp(mode,'generate')
%  This just adds this element to the element data structure. 
  elnum=c;
  element(elnum).nodes=b(1);
%  element(elnum).properties=b(2);
end
if strcmp(mode,'make')
% There is nothing to make with a surface node element. It's
% existence is everything it does.
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end
