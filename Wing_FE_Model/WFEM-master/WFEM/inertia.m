function out=inertia(mode,b,c,d,e)

% INERTIA adds point masses and inertias.
% The inertia properties can be entered in the following forms and
% are a "material type". 
% For a sphere:
% m I
% For an irregular object:  
% m Ixx Ixy Ixz Iyy Iyz Izz
% 
% To be added at a later date:
% For an irregular object not oriented with the global x-y-z
% coordinates:
% m Ixx Ixy Ixz Iyy Iyz Izz l1x ly1 lz1 l2x l2y l2z 
% where the local x coordinate points in the direction of the
% vector [l1x l1y l1z], and the local y coordinate points in the
% direction of the vector [l2x l2y l2z]
% 
% Defining inertia element in wfem input file:
% node materialnumber
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
mode;
if strcmp(mode,'numofnodes')
  out=1;
elseif strcmp(mode,'nastrandump')% No properties needed |strcmp(mode,'nastranpropertydump')
  elnum=b;
  fido=c;
  node=element(elnum).nodes;
  iprops=elprops(element(elnum).properties).a;% element(elnum).properties 
                                              % stores the
                                              % properties number
                                              % of the current
                                              % elnum. elprops
                                              % contains this data.
  iprops;
  m=iprops(1);
  if length(iprops)==1
    I=zeros(3);
  elseif length(iprops)==2
    I=eye(3)*iprops(2);
  elseif length(iprops)==7
    Ixx=iprops(2);
    Ixy=iprops(3);
    Ixz=iprops(4);
    Iyy=iprops(5);
    Iyz=iprops(6);
    Izz=iprops(7);
    I=[Ixx Ixy Ixz;...
       Ixy Iyy Iyz;...
       Ixz Iyz Izz];
  else
    disp('Poorly defined material property for inertia element')
    I=zeros(3);
  end
  Ixx=I(1,1);
  Ixy=I(1,2);
  Ixz=I(1,3);
  Iyy=I(2,2);
  Iyz=I(2,3);
  Izz=I(3,3);
  fprintf(fido,'CONM2, %g, %g, 0 , %8.4E, 0.0, 0.0, 0.0 \n , %8.4E,%8.4E,%8.4E,%8.4E,%8.4E,%8.4E \n',...
	  [elnum, node, m, Ixx, Ixy, Iyy, Ixz, Iyz, Izz])
  
elseif strcmp(mode,'generate')
  elnum=c;
  element(elnum).nodes=b(1);
  element(elnum).properties=b(2);
elseif strcmp(mode,'make')
  elnum=b;%element(elnum)
  node=element(elnum).nodes; 
  iprops=elprops(element(elnum).properties).a;% element(elnum).properties 
                                              % stores the
                                              % properties number
                                              % of the current
                                              % elnum. elprops
                                              % contains this data.
  iprops;
  m=iprops(1);
  if length(iprops)==1
    I=zeros(3);
  elseif length(iprops)==2
    I=eye(3)*iprops(2);
  elseif length(iprops)==7
    Ixx=iprops(2);
    Ixy=iprops(3);
    Ixz=iprops(4);
    Iyy=iprops(5);
    Iyz=iprops(6);
    Izz=iprops(7);
    I=[Ixx Ixy Ixz;...
       Ixy Iyy Iyz;...
       Ixz Iyz Izz];
  else
    disp('Poorly defined material property for inertia element')
    I=zeros(3);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  indices=[(1:6)+(node-1)*6];
  %indices
  %m
  %I
  %M,pause
  M(indices,indices)=M(indices,indices)+[m*eye(3) zeros(3);zeros(3) ...
		    I];
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end
