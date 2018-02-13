function out=sensor(mode,b,c,d,e)

% SENSOR denotes a simple. No calibration/errors are possible. All
% three DOFs of the node are presumed to be sensed.
%  
% SENSOR properties are simple:
% 1: displacement
% 2: velocity
% 3L acceleration
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables (global):
% -------------------
% K       :    Global stiffness matrix
% Ks      :    Global stiffness buckling matrix
% M       :    Global mass matrix
% nodes   :    [x y z] nodal locations
% Cd      :    Displacement Measurement
% Cv      :    Velocity Measurement
% Ca      :    Acceleration Measurement
% Xd      :    Original Positions  
% nout    :    number of outputs defined currently (of any type)
  
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
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Joseph C. Slater, 11/28/2005.
% joseph.slater@wright.edu
out=0;
if strcmp(mode,'numofnodes')
  out=1;
end
if strcmp(mode,'generate')
%  This just adds this element to the element data structure. 
  elnum=c;
  
  element(elnum).nodes=b(1) ;
  element(elnum).properties=b(2);
  element(elnum).lineno=curlineno;
  %How should I draw the darn thing. Likely a surface. Right. Using
  %the draw call, or not?
  
end

if strcmp(mode,'make')
  elnum=b;
  elnodes=element(elnum).nodes;
  senstype=elprops(element(elnum).properties).a;
  indices=[elnodes*6+(-5:0) ];

  % There is nothing to make with a surface node element. Its
  % existence is everything it does.
  
  if senstype==1
    scd=size(Cd,1);
    Cd(nout+(1:6),elnodes*6+(-5:0))=eye(6);
  elseif senstype==2
    scv=size(Cd,1);
    Cv(nout+(1:6),elnodes*6+(-5:0))=eye(6);
  elseif senstype==3
    sca=size(Cd,1);
    Cd(nout+(1:6),elnodes*6+(-5:0))=eye(6);
  else
    disp(['Warning: Unknown sensor number ' num2str(senstype) '.'])
    pause
  end
  
  nout=nout+6;
  
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end
