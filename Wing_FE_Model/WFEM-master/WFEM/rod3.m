function out=rod3(mode,b,c,d,e)
  
% ROD3 does as listed below. It is a simple 2-nodes rod element in 3-D
% The rod properties (bprops) are in the order
% bprops=[E rho A]
% All formats other than this are reduced to these properties by
% averaging or discarding.
%
% Defining beam3 element properties in wfem input file:
% element properties
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 
%   E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2 
%   E G rho A J Izz Iyy 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...  
%       mx2 my2 mz2 mrx2 mry2 mrz2 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...  
%       sx2 sy2 sz2 srx2 sry2 srz2 distype 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...  
%       mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 distype 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
%       mx2 my2 mz2 mrx2 mry2 mrz3 mx3 my3 mz3 mrx3 mry3 mrz3 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
%       mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...  
%       mx3 my3 mz3 mrx3 mry3 mrz3 sx3 sy3 sz3 srx3 sry3 srz3 distype
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
%       mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...  
%       lx2 ly2 lz2 lrx2 lry2 lrz2 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
%       mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...  
%       lx2 ly2 lz2 lrx2 lry2 lrz2 ...  
%       mx3 my3 mz3 mrx3 mry3 mrz3 sx3 sy3 sz3 srx3 sry3 srz3 ...
%       lx3 ly3 lz3 lrx3 lry3 lrz3 
%   E rho A 
%   E rho A m
%   E rho A m s
%   E rho A m s distype
%   E rho A m l distype
%   E rho A m s l 
% Property formats other than E rho A (plus additional strain
% options) add no capabilities, but are available
% only for substituting rod3 elements for beam3 elements. m, s, l
% stand for mean, standard deviation, or limit. A distype of zero
% means Gaussian. A distype of -1 means uniform. 
%
% Defining rod3 element in wfem input file:
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
% Copyright Joseph C. Slater, 9/1/2003.
% joseph.slater@wright.edu
out=0;
if strcmp(mode,'numofnodes')
  % return the number of modes in the variable 'out'. We are
  % accommodating beam3 format. 
  if length(b)==5, 
    out=3;
  else
    out=2;
  end
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
	  'A rod element needs 2 nodes and a material number. \n'...
	  'Assuming material number 1.']))
    beep
  else
    propnumber=b(length(b));
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
  G=bprops(2);
  rho=bprops(3);
  ss=0;mm=0; %standard deviation and mean of initial deformation defaults
  if length(bprops)<=6
    rho=bprops(2);
    A=bprops(3);
    if length(bprops)>3
      mm=bprops(4);
    end
    if length(bprops)==5
      ss=bprops(5);
      distype=0;
    end
    if length(bprops)==6
      if bprops(6)<0
	distype=bprops(6);
	ss=bprops(5);
	ll=bprops(5);
      else
	ss=bprops(5);
	ll=bprops(6);
	distype=-1;
      end
    end
    
    %  length(bprops),pause
  elseif length(bprops)==15
    A=(bprops(4)+bprops(5)+bprops(6))/3;
  elseif length(bprops)==11
    A=(bprops(4)+bprops(5))/2;
  elseif length(bprops)==7
    A=bprops(4);
  elseif length(bprops)==21
    A=(bprops(4)+bprops(5)+bprops(6))/3;
    mm=bprops(16);
    ss=0;
    distype=-1;
  elseif length(bprops)==27
    A=(bprops(4)+bprops(5)+bprops(6))/3;
    mm=bprops(16);
    ss=0;
  elseif length(bprops)==28
    A=(bprops(4)+bprops(5)+bprops(6))/3;
    mm=bprops(16);
    ss=bprops(22);
    distype=-abs(bprops(28));
  elseif length(bprops)==40
    A=(bprops(4)+bprops(5)+bprops(6))/3;
    mm=bprops(16);
    ss=bprops(22);
    distype=-abs(bprops(40));
  elseif length(bprops)==51|length(bprops)==33
    A=(bprops(4)+bprops(5)+bprops(6))/3;
    mm=bprops(16);
    ss=bprops(22);
    ll=bprops(28);
  else
    warndlg(['The number of material properties set for ' ...
             'this element (' num2str(length(bprops)) ') isn''t ' ...
             'appropriate for a beam3 element. '      ...
             'Please refer to the manual.'],...
             'Bad element property definition.','modal');
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
  if sqrt(A)/ellength<.002
    disp(['Element ' num2str(elnum) ' is exceptionally thin.'...
	     'Please be aware.'])
  end
  
  
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
