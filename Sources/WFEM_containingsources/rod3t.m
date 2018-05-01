function out=rod3t(mode,b,c,d,e)
  
% ROD3T is a simple 2-nodes rod element in 3-D with both thermal
% strain and initial strain (length error) capability. 
% The rod properties (bprops) are in the order
% bprops=[E rho A]
% All formats other than this are reduced to these properties by
% averaging or discarding.
%
% Defining rod3t element properties in wfem input file:
% element properties
%   E rho A m s distype alpha alim thermalinput
%   E rho A m l distype alpha alim thermalinput
%   E rho A m s l       alpha alim thermalinput
% m, s, l stand for mean, standard deviation, or limit. A distype of zero
% means Gaussian. A distype of -1 means uniform. 
% alpha is thermal coefficient of expansion. alim is the upper
% uniform variation bound (delta alpha max). thermalinput defines
% which temperature input effects this strut.
%
% Defining rod3t element in wfem input file:
% node1 node2 materialnumber
%
% See wfem.m for more explanation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables (global):
% -------------------
% K       :    Global stiffness matrix
% Ks      :    Global stiffness buckling matrix
% M       :    Global mass matrix
% nodes   :    [x y z] nodal locations
global ismatnewer
global K
global Ks
global M
global nodes % Node locations
global elprops
global element
global points
global Fepsn % Initial strain "forces". 
global lines
global restart
global reload
global curlineno
global Bthermal %thermal input matrix
%
% Variables (local):
% ------------------
% bnodes  :    node/point numbers for actual beam nodes 1-2-3 and point
% k       :    stiffness matrix in local coordiates
% kg      :    stiffness matrix rotated into global coordinates
% m       :    mass matrix in local coordiates
% mg      :    mass matrix rotated into global coordinates
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Joseph C. Slater, 9/1/2005.
% joseph.slater@wright.edu
out=0;
if strcmp(mode,'numofnodes')
  % return the number of modes in the variable 'out'. 
  out=2;
elseif strcmp(mode,'generate')
  elnum=c;
  % The second argument (b) is the line of material properties read
  % from the input file. We're going to start filing them  away in
  % more meaningful variables. 
  
  element(elnum).nodes=[b(1) b(2)];% I don't care how long the
				   % list is, for a beam3 or a
				   % rod, the first two items are
				   % the end nodes.
  if length(b)==2
    propnumber=1;
    disp(sprintf(['Insufficient definition of rod element. \n' ...
	  'A rod3t element needs 2 nodes and a material number. \n'...
	  'Assuming material number 1.']))
    beep
  else
    propnumber=b(3);
  end
  
  element(elnum).properties=propnumber;
  element(elnum).lineno=curlineno; % We record this so we can look
                                   % it up later for error messages.
end

% Here we figure out what the rod properties mean. If you need
% them in a call-mode, that mode should be in the if on the next line.
if strcmp(mode,'make')|strcmp(mode,'istrainforces')
  elnum=b;element(elnum);
  bnodes=[element(elnum).nodes]; % This is shorter to type, and
                                 % won't conflict with the nodes
                                 % variable already in use.
%The following comment section left in for the sake of remembering
%some things. It should be removed later.
                     %element(elnum).point];% The point is
                                                     % referred to
                                                     % as node 4
                                                     % below,
                                                     % although it
                                                     % actually
                                                     % calls the
                                                     % array points
                                                     % to get it's
                                                     % location.
% elnum
% disp('property number')
% element(elnum).properties
% elprops

  bprops=elprops(element(elnum).properties).a;% element(elnum).properties 
                                              % stores the
                                              % properties number
                                              % of the current
                                              % elnum. elprops
                                              % contains this data.

  E=bprops(1);
  rho=bprops(2);
  A=bprops(3);
  ss=0;mm=0; %standard deviation and mean of initial deformation
             %defaults
  alpha=bprops(7);
  alim=bprops(8);
  mm=bprops(4);
  
  
  thermalinput=bprops(9);
  if bprops(6)==-1
    %uniform distribution
    ll=bprops(5);distype=-1;
  elseif bprops(6)==0
    %gaussian distribution
    ss=bprops(5);distype=0;
  else
    %Gaussian truncated distribution
    ss=bprops(5);
    ll=bprops(6);
  end
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rod properties (bprops) are in the order
% bprops=[E rho A]

if strcmp(mode,'make')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Define beam node locations for easy later referencing
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  x1=nodes(bnodes(1),1);
  y1=nodes(bnodes(1),2);
  z1=nodes(bnodes(1),3);
  x2=nodes(bnodes(2),1);
  y2=nodes(bnodes(2),2);
  z2=nodes(bnodes(2),3);

  % This rod code is 'hard coded', including rotations, for
  % speed. To see a more general example, please see beam3.m
  
  % Modified notation from Rao, "The Finite Element Method in Engineering",
  % 3rd Edition.
  l=x2-x1;
  m=y2-y1;
  n=z2-z1;
  ellength=norm([l m n]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now we create the thermal disturbance columns for this element.
  
  %Apply uniform distribution
  alpha=alpha+(rand(1)-.5)*alim*2;
  
  %Create thermal force vector
  tdirvect=[-l;-m;-n;l;m;n]/ellength;
  tvect=E*A*alpha*tdirvect;
  
  

% $$$   if sqrt(A)/ellength<.002
% $$$     disp(['Element ' num2str(elnum) ' is exceptionally thin.'...
% $$$ 	     'Please be aware.'])
% $$$   end
  
  
  ellengthsq=ellength^2;
  premul=E*A/ellengthsq/ellength;
  ll=premul*l^2;
  lm=premul*l*m;
  ln=premul*l*n;
  mm=premul*m^2;
  mn=premul*m*n;
  nn=premul*n^2;
  k=[ ll  lm  ln -ll -lm -ln;
      lm  mm  mn -lm -mm -mn;
      ln  mn  nn -ln -mn -nn;
     -ll -lm -ln  ll  lm  ln;
     -lm -mm -mn  lm  mm  mn;
     -ln -mn -nn  ln  mn  nn];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % Derivation of Mass matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  two=rho*A*ellength/3;
  one=two/2;
  m=[two 0   0   one 0   0  ;
     0   two 0   0   one 0  ;
     0   0   two 0   0   one;
     one 0   0   two 0   0  ;
     0   one 0   0   two 0  ;
     0   0   one 0   0   two];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Store elemental matrices in case we need them later
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  element(elnum).m=m;
  element(elnum).k=k;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  bn1=bnodes(1);bn2=bnodes(2);
  indices=[bn1*6+(-5:-3) bn2*6+(-5:-3)] ;%Don't forget that the
                                       %global system has 6 DOF at
                                       %each node, while our
                                       %element only has 3.

  K(indices,indices)=K(indices,indices)+k;
  M(indices,indices)=M(indices,indices)+m;

  if size(Bthermal,2)<thermalinput
    Bthermal(size(Bthermal,1),thermalinput)=0;
  end
  
  %size(Bthermal)
  %thermalinput
  %pause
  Bthermal(indices,thermalinput)=Bthermal(indices,thermalinput)+tvect;
  
  
  % At this point we also know how to draw the element (what lines
  % and surfaces exist). For the beam3 element, 2 lines are
  % appropriate.
  numlines=size(lines,1);
  lines(numlines+1,:)=[bn1 bn2];
  
%diag(M)
elseif strcmp(mode,'istrainforces')
  % We need to have the stiffness matrix and the coordinate roation matrix.
  bn1=bnodes(1);bn2=bnodes(2);
  x1=nodes(bn1,1);
  y1=nodes(bn1,2);
  z1=nodes(bn1,3);
  x2=nodes(bn2,1);
  y2=nodes(bn2,2);
  z2=nodes(bn2,3);
  l=x2-x1;
  m=y2-y1;
  n=z2-z1;
  ellengthsq=l^2+m^2+n^2;
  ellength=sqrt(ellengthsq);

  k=element(elnum).k;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Add initial deflections
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if exist('distype')==0
    distype=0;
  end
  
  if distype==0% Gaussian distribution
	       %disp('need to check for limits, then apply them/make them if
	       %they don't exist')
	       %normal distribution
    if ~exist('ll')
      ll=inf;
    end
    dev=randn(1)*ss;
    while dev>ll
      dev=randn(1)*ss;
    end
    dx=mm+dev;
  elseif distype==-1
    %uniform distribution
    dx=mm+(rand(1)-.5)*ll*2;
  else
    disp(['Invalid distribution type ' num2str(distype) '.'])
  end

  % We're doing everything in global coordinates to save time. 
  Fistrain=k*[0;0;0;
	      [l;m;n]*dx/ellength];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  indices=[bn1*6+(-2:0) bn2*6+(-2:0)] ;
  
  if length(Fepsn)<max(indices);
    Fepsn(max(indices))=0;
  end
  Fepsn(indices)=Fepsn(indices)+Fistrain;

  
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end
