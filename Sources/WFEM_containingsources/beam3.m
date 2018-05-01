function out=beam3(mode,b,c,d,e)
  
% BEAM3 does as listed below. It is an Euler-Bernoulli
% beam/rod/torsion model. 
% Beam properties (bprops) are in the order
% bprops=[E G rho A1 A2 A3 J1 J2 J3 Ixx1 Ixx2 Ixx3 Iyy1 Iyy2 Iyy3]
% For a linear interpolation they are
% bprops=[E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2]
% Note that the linear interpolation a user shortcut
% and results in no less computational effort.
% Third node is in the middle.
% Fourth "node" defines the beam y plane and is actually from the
% points array.
%
%
% Defining beam element properties in wfem input file:
% element properties
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 
%   E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2 
%   E G rho A J Izz Iyy 
%   E G rho A J Izz Iyy sx2 sy2 sz2 srx2 sry2 srz2 distype 
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
% If the second format is used, a linear interpolation of the
% properties is presumed. If the third format is used, constant
% properties are assumed. In subsequent lines, initial deflections
% can be prescribed, deterministically, or stochastically
% (random). The character \command{m} stands for \emph{mean} value,
% and \command{r} stands for \emph{rotation} value. If a stochastic
% form is used, a distribution type must be prescribes. A
% \command{distype} of \command{0} means normal distribution with
% standard deviation of \command{s} values as prescribe. A
% \command{distype} of \command{-1} means uniform distribution with
% bounds relative to the mean set by the \command{s} values. A
% truncated Gaussian distribution is demonstrated in the last
% line. There the \command{l} values are limits relative to the
% mean (just as the \command{s} values for a uniform
% distribution). Note that the dots illustrate continuation of the
% same line. Continuations of a line using the M\textsc{atlab}
% notation \ldots can only be used in defining element
% properties. Torsional rigidity, $J$, must be less than or equal
% to $Iyy+Izz$ at any given cross section.  
%
% Defining beam3 element in wfem input file:
%   node1 node2 node3 pointnumber materialnumber 
%   node1 node2 pointnumber materialnumber 
%   node1 node2 materialnumber
%
% node2 is automatically created if it is missing.
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
global DoverL
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
% Copyright Joseph C. Slater, 7/26/2002.
% joseph.slater@wright.edu
out=0;
if strcmp(mode,'numofnodes')
  if length(b)==5, 
    out=3;
  else
    out=2;
  end
end
if strcmp(mode,'generate')
  elnum=c;b;
  if length(b)==3|length(b)==4
    %nodes,pause
    n3=size(nodes,1)+1;
    element(elnum).nodes=[b(1) b(2) n3];% recall node 3 is in the
                                        % middle.
    %exist('d')
    %if exist('d')==0
      nodes(n3,1:3)=(nodes(b(1),:)+nodes(b(2),:))/2;% We
						    % need
						    % to
						    % add
						    % this
						    % as a node, but
						    % not if we are
						    % just generating
						    % a nastran
                                                    % input
                                                    % deck. OK,
                                                    % Greg changed
                                                    % his mind,
                                                    % we'll make it
    %end
    
    
    if length(b)==3
      element(elnum).properties=b(3);
      element(elnum).point=1;
    else
      element(elnum).point=b(3);
      element(elnum).properties=b(4);
    end
  elseif length(b)==5
    element(elnum).nodes=b(1:3);
    if norm((nodes(b(1),:)+nodes(b(2),:))/2-nodes(b(3),:))/ ...
	  norm(nodes(b(1),:)-nodes(b(2),:))>.001
      disp(['WARNING: Node ' num2str(b(3)) ...
	    ' is not in the middle of element ' ...
	    num2str(elnum) ' on line ' num2str(curlineno) '.'])
    end
    element(elnum).properties=b(5);
    element(elnum).point=b(4);
  else 
    warndlg(['Element ' num2str(elnum) ' on line ' ...
	     num2str(element(elnum).lineno) ' entered incorrectly.'], ...
            ['Malformed Element'],'modal')
    return
  end
  element(elnum).lineno=curlineno;
  element(elnum);
  if exist('d')==1% If this is going to be making nastran elements,
                  % we need to split it into 2 2-noded elements. 
      elnum2=length(element)+1;
      element(elnum2)=element(elnum);
      element(elnum2).nodes(1)=element(elnum2).nodes(3);
      element(elnum).nodes(2)=element(elnum).nodes(3);
  end
  
end

if strcmp(mode,'nastrandump')|strcmp(mode,'nastranpropertydump')
  elnum=b;
  fido=c;
  node1=element(elnum).nodes(1);
  node2=element(elnum).nodes(2);
  propnum=element(elnum).properties;
  point=element(elnum).point;
  if strcmp(mode,'nastrandump')
    fprintf(fido,'CBAR,%g, %g, %g, %g, %8.4E ,%8.4E,%8.4E \n',[elnum, propnum, ...
		    node1, node2, 1,0,0,]);
    %1, 0, 0 is the point required by FEMLAB. It's pretty useless
    %still. 
% $$$     fprintf(fido,['CBAR,' num2str(elnum) ', ' num2str(propnum) ', ' ...
% $$$ 		  num2str(node1) ', ' num2str(node2) ', ' num2str(point) ...
% $$$ 		  '\n']);
  else
    %We should move inside the next loop and average. Also note
    %that the use of the point hasn't been checked to match that of
    %NASTRAN. Symmetry is allowable only. 
    bprops=elprops(element(elnum).properties).a;
    E=bprops(1);
    G=bprops(2);
    rho=bprops(3);
    A=bprops(4);
    Iyy=bprops(7);
    Izz=bprops(6);
    J=bprops(5);
    
    fprintf(fido,'MAT1,%g,%8.4E,%8.4E,,%8.4E,0.001,0.001,0.001 \n',[propnum,E,G,rho]);
% $$$     fprintf(fido,['MAT1,' num2str(propnum) ', '  num2str(E) ...
% $$$ 		  ', ' num2str(G) ', ,' num2str(rho) ', ' ...
% $$$ 		  num2str(.001) ', ' num2str(.001) ', ' num2str(.001) '\n']);
    fprintf(fido,'PBAR, %g, %g, %8.4E, %8.4E, %8.4E, %8.4E, 0\n',[propnum, ...
		    propnum, A, Iyy, Izz, J]);
% $$$     fprintf(fido,['PBAR,' num2str(propnum) ', '  num2str(propnum) ', ' num2str(A) ...
% $$$ 		  ', ' num2str(Iyy) ', ' num2str(Izz) ',' num2str(J) ...
% $$$ 		  ',0\n' ]);
  end
  
end





% Here we figure out what the beam properties mean. If you need
% them in a mode, that mode should be in the if on the next line.
if strcmp(mode,'make')|strcmp(mode,'istrainforces')
  elnum=b;element(elnum);
  bnodes=[element(elnum).nodes element(elnum).point];% The point is
                                                     % referred to
                                                     % as node 4
                                                     % below,
                                                     % although it
                                                     % actually
                                                     % calls the
                                                     % array points
                                                     % to get it's location.
  bprops=elprops(element(elnum).properties).a;% element(elnum).properties 
                                              % stores the
                                              % properties number
                                              % of the current
                                              % elnum. elprops
                                              % contains this data.
  
  E=bprops(1);
  G=bprops(2);
  
  rho=bprops(3);
%  length(bprops),pause
  if length(bprops)==15
    A1=bprops(4);
    A2=bprops(5);
    A3=bprops(6);
    J1=bprops(7);
    J2=bprops(8);
    J3=bprops(9);
    Izz1=bprops(10);
    Izz2=bprops(11);
    Izz3=bprops(12);
    Iyy1=bprops(13);
    Iyy2=bprops(14);
    Iyy3=bprops(15);
  elseif length(bprops)==11
    A1=bprops(4);
    A2=bprops(5);
    A3=(A1+A2)/2;
    J1=bprops(6);
    J2=bprops(7);
    J3=(J1+J2)/2;
    Izz1=bprops(8);
    Izz2=bprops(9);
    Izz3=(Izz1+Izz2)/2;
    Iyy1=bprops(10);
    Iyy2=bprops(11);  
    Iyy3=(Iyy1+Iyy2)/2;
  elseif length(bprops)==7
    A1=bprops(4);
    A2=A1;A3=A1;
    J1=bprops(5);
    J2=J1;J3=J1;
    Izz1=bprops(6);
    Izz2=Izz1;Izz3=Izz1;
    Iyy1=bprops(7);
    Iyy2=Iyy1;Iyy3=Iyy1;
  elseif length(bprops)==14
    A1=bprops(4);
    A2=A1;A3=A1;
    J1=bprops(5);
    J2=J1;J3=J1;
    Izz1=bprops(6);
    Izz2=Izz1;Izz3=Izz1;
    Iyy1=bprops(7);
    Iyy2=Iyy1;Iyy3=Iyy1;
    distype=bprops(14);
    if distype==0
      sdx2=bprops(8);
      sdy2=bprops(9);
      sdz2=bprops(10);
      sdrx2=bprops(11);
      sdry2=bprops(12);
      sdrz2=bprops(13);
      mdx2=0;
      mdy2=0;
      mdz2=0;
      mdrx2=0;
      mdry2=0;
      mdrz2=0;
    else
      sdx2=bprops(8);
      sdy2=bprops(9);
      sdz2=bprops(10);
      sdrx2=bprops(11);
      sdry2=bprops(12);
      sdrz2=bprops(13);
      mdx2=0;
      mdy2=0;
      mdz2=0;
      mdrx2=0;
      mdry2=0;
      mdrz2=0;
    end
    
  elseif length(bprops)==21
    A1=bprops(4);
    A2=bprops(5);
    A3=bprops(6);
    J1=bprops(7);
    J2=bprops(8);
    J3=bprops(9);
    Izz1=bprops(10);
    Izz2=bprops(11);
    Izz3=bprops(12);
    Iyy1=bprops(13);
    Iyy2=bprops(14);
    Iyy3=bprops(15);
    mdx2=bprops(16);
    mdy2=bprops(17);
    mdz2=bprops(18);
    mdrx2=bprops(19);
    mdry2=bprops(20);
    mdrz2=bprops(21);
    sdx2=0;
    sdy2=0;
    sdz2=0;
    sdrx2=0;
    sdry2=0;
    sdrz2=0;
    disp('Uniform node 2 deflections.')
  elseif length(bprops)==27
    A1=bprops(4);
    A2=bprops(5);
    A3=bprops(6);
    J1=bprops(7);
    J2=bprops(8);
    J3=bprops(9);
    Izz1=bprops(10);
    Izz2=bprops(11);
    Izz3=bprops(12);
    Iyy1=bprops(13);
    Iyy2=bprops(14);
    Iyy3=bprops(15);
    mdx2=bprops(16);
    mdy2=bprops(17);
    mdz2=bprops(18);
    mdrx2=bprops(19);
    mdry2=bprops(20);
    mdrz2=bprops(21);
    mdx3=bprops(22);
    mdy3=bprops(23);
    mdz3=bprops(24);
    mdrx3=bprops(25);
    mdry3=bprops(26);
    mdrz3=bprops(27);
    sdx2=0;
    sdy2=0;
    sdz2=0;
    sdrx2=0;
    sdry2=0;
    sdrz2=0;
    sdx3=0;
    sdy3=0;
    sdz3=0;
    sdrx3=0;
    sdry3=0;
    sdrz3=0;
  elseif length(bprops)==28
    A1=bprops(4);
    A2=bprops(5);
    A3=bprops(6);
    J1=bprops(7);
    J2=bprops(8);
    J3=bprops(9);
    Izz1=bprops(10);
    Izz2=bprops(11);
    Izz3=bprops(12);
    Iyy1=bprops(13);
    Iyy2=bprops(14);
    Iyy3=bprops(15);
    mdx2=bprops(16);
    mdy2=bprops(17);
    mdz2=bprops(18);
    mdrx2=bprops(19);
    mdry2=bprops(20);
    mdrz2=bprops(21);
    sdx2=bprops(22);
    sdy2=bprops(23);
    sdz2=bprops(24);
    sdrx2=bprops(25);
    sdry2=bprops(26);
    sdrz2=bprops(27);
    distype=bprops(28);
  elseif length(bprops)==40
    A1=bprops(4);
    A2=bprops(5);
    A3=bprops(6);
    J1=bprops(7);
    J2=bprops(8);
    J3=bprops(9);
    Izz1=bprops(10);
    Izz2=bprops(11);
    Izz3=bprops(12);
    Iyy1=bprops(13);
    Iyy2=bprops(14);
    Iyy3=bprops(15);
    mdx2=bprops(16);
    mdy2=bprops(17);
    mdz2=bprops(18);
    mdrx2=bprops(19);
    mdry2=bprops(20);
    mdrz2=bprops(21);
    sdx2=bprops(22);
    sdy2=bprops(23);
    sdz2=bprops(24);
    sdrx2=bprops(25);
    sdry2=bprops(26);
    sdrz2=bprops(27);
    mdx3=bprops(28);
    mdy3=bprops(29);
    mdz3=bprops(30);
    mdrx3=bprops(31);
    mdry3=bprops(32);
    mdrz3=bprops(33);
    sdx3=bprops(34);
    sdy3=bprops(35);
    sdz3=bprops(36);
    sdrx3=bprops(37);
    sdry3=bprops(38);
    sdrz3=bprops(39);
    distype=bprops(40);
  elseif length(bprops)==51|length(bprops)==33
    A1=bprops(4);
    A2=bprops(5);
    A3=bprops(6);
    J1=bprops(7);
    J2=bprops(8);
    J3=bprops(9);
    Izz1=bprops(10);
    Izz2=bprops(11);
    Izz3=bprops(12);
    Iyy1=bprops(13);
    Iyy2=bprops(14);
    Iyy3=bprops(15);
    mdx2=bprops(16);
    mdy2=bprops(17);
    mdz2=bprops(18);
    mdrx2=bprops(19);
    mdry2=bprops(20);
    mdrz2=bprops(21);
    sdx2=bprops(22);
    sdy2=bprops(23);
    sdz2=bprops(24);
    sdrx2=bprops(25);
    sdry2=bprops(26);
    sdrz2=bprops(27);
    ldx2=bprops(28);
    ldy2=bprops(29);
    ldz2=bprops(30);
    ldrx2=bprops(31);
    ldry2=bprops(32);
    ldrz2=bprops(33);
    if length(bprops)==51
      % More complete definition of the initial deflections due to
      % initial strain. 
      mdx3=bprops(34);
      mdy3=bprops(35);
      mdz3=bprops(36);
      mdrx3=bprops(37);
      mdry3=bprops(38);
      mdrz3=bprops(39);
      sdx3=bprops(40);
      sdy3=bprops(41);
      sdz3=bprops(42);
      sdrx3=bprops(43);
      sdry3=bprops(44);
      sdrz3=bprops(45);
      ldx3=bprops(46);
      ldy3=bprops(47);
      ldz3=bprops(48);
      ldrx3=bprops(49);
      ldry3=bprops(50);
      ldrz3=bprops(51);
    end
  else
    warndlg(['The number of material properties set for ' ...
             'this element (' num2str(length(bprops)) ') isn''t ' ...
             'appropriate for a beam3 element. '      ...
             'Please refer to the manual.'],...
             'Bad element property definition.','modal');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beam properties (bprops) are in the order
% bprops=[E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3]
% For a linear beam they are
% bprops=[E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2]

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
  x3=nodes(bnodes(3),1);
  y3=nodes(bnodes(3),2);
  z3=nodes(bnodes(3),3);
  x4=points(bnodes(4),1);
  y4=points(bnodes(4),2);
  z4=points(bnodes(4),3);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Shape functions for higher order beam. 
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Shape functions in matrix polynomial form (polyval style) for bending
  bn1 =  [  0.750  -0.500  -1.250   1.000   0.000   0.000];
  bn1d =  [3.75000  -2.00000  -3.75000   2.00000   0.00000];
  bn1dd =  [   15.00   -6.00   -7.50    2.00];
  bn2 =  [ 0.250  -0.250  -0.250   0.250   0.000   0.000];
  bn2d =  [1.25000  -1.00000  -0.75000   0.50000   0.00000];
  bn2dd =  [   5.000  -3.000  -1.500   0.500];
  bn3 =  [-0.750  -0.500   1.250   1.000   0.000   0.000];
  bn3d = [-3.75000  -2.00000   3.75000   2.00000   0.00000];
  bn3dd =[  -15.00   -6.00    7.50    2.00];
  bn4 =  [ 0.250   0.250  -0.250  -0.250   0.000   0.000];
  bn4d = [ 1.25000   1.00000  -0.75000  -0.50000   0.00000];
  bn4dd =[   5.000   3.000  -1.500  -0.500];
  bn5 =  [ 0.000   1.000  -0.000  -2.000   0.000   1.000];
  bn5d = [0.00000   4.00000  -0.00000  -4.00000   0.00000];
  bn5dd =[   0    1.20e+01   0   -4.00e+00];
  bn6 =  [ 1.000   0.000  -2.000  -0.000   1.000   0.000];
  bn6d = [ 5.0000    0.0   -6.0000   0.0    1.0000];
  bn6dd =[    20.0    0   -12.0  0];
  
  % Shape functions in matrix polynomial form (polyval style) for 
  % torsion/rod
  rn1=[0.5 -.5 0];
  rn1d=[1 -0.5];
  rn2=[.5 .5 0];
  rn2d=[1 0.5];
  rn3=[-1 0 1];
  rn3d=[-2 0];
  numbeamgauss=5; % Number of Gauss points for integration of beam element
  [bgpts,bgpw]=gauss(numbeamgauss);
  kb1=zeros(6,6);
  kb2=kb1;
  l=norm([x2 y2 z2]-[x1 y1 z1]);
  propertynum=num2str(element(elnum).properties);
  % Allowable aspect ratio. I recommend D/l=.1
  if isempty(DoverL)==1
    DoverL=.1;
  end
  
  if sqrt(A1*4/pi)/l>DoverL|sqrt(A2*4/pi)/l>DoverL|sqrt(A3*4/pi)/l>DoverL
    warndlg({['Dimensions of element ' num2str(elnum) ' using properties '...
	      propertynum ' are more suitable for a Timoshenko beam.'];...
	     'radius divided by length is too large'},...
	    'Improper application of element.','replace')
  end
  if (Izz1+Iyy1)<(1/2.1*A1^2/pi)|(Izz2+Iyy2)<(1/2.1*A2^2/pi)|(Izz3+Iyy3)<(1/2.1*A3^2/pi)
    %2.0 would be exact for a circle
    warndlg({['Iyy+Izz for properties number' propertynum ' can''t be as '...
	      'low as have been given.'];...
	     'Nonphysical properties.'},['Impossible cross sectional' ...
		    ' properties'],'replace')
  end
  slenderness=min([sqrt((Izz1+Iyy1)/A1) sqrt((Izz2+Iyy2)/A2) ...
		   sqrt((Izz3+Iyy3)/A3)  ])/l;
  if slenderness<.002
    disp([num2str(elnum) ['is a rediculously thin element. Please' ...
		    ' check numbers.']])
  end
  
  Jac=l/2;% Beam Jacobian. valid only if node three is in the
          % middle of the beam
          % Local Bending in x-y plane
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1dd,bgpts(i))/Jac^2;
             polyval(bn2dd,bgpts(i))/Jac;
             polyval(bn3dd,bgpts(i))/Jac^2;
             polyval(bn4dd,bgpts(i))/Jac;
             polyval(bn5dd,bgpts(i))/Jac^2;
             polyval(bn6dd,bgpts(i))/Jac];
    Izz=polyval(rn1*Izz1+rn2*Izz2+rn3*Izz3,bgpts(i));%these should
                                                     %be called Izz
    kb1=kb1+bgpw(i)*beamsfs*beamsfs'*Izz*E*Jac;
  end
  % Local Bending in x-z plane
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1dd,bgpts(i))/Jac^2;
             -polyval(bn2dd,bgpts(i))/Jac;
             polyval(bn3dd,bgpts(i))/Jac^2;
             -polyval(bn4dd,bgpts(i))/Jac;
             polyval(bn5dd,bgpts(i))/Jac^2;
             -polyval(bn6dd,bgpts(i))/Jac];
    Iyy=polyval(rn1*Iyy1+rn2*Iyy2+rn3*Iyy3,bgpts(i));
    kb2=kb2+bgpw(i)*beamsfs*beamsfs'*Iyy*E*Jac;
  end
  
  % Local Extension in x, torsion about x
  numrodgauss=3;
  [rgpts,rgpw]=gauss(numrodgauss);
  krod=zeros(3,3);
  ktor=zeros(3,3);
  for i=1:numrodgauss
    rodsfs=[polyval(rn1d,rgpts(i))/Jac;
            polyval(rn2d,rgpts(i))/Jac;
            polyval(rn3d,rgpts(i))/Jac];
    if (J1>(Iyy1+Izz1))|(J2>(Iyy2+Izz2))|(J3>(Iyy3+Izz3))
      if (J1>(Iyy1+Izz1))
	disp('WARNING: J1 must be <= Iyy1+Izz1')
      end
      if (J2>(Iyy2+Izz2))
	disp('WARNING: J2 must be <= Iyy2+Izz2')
      end
      if (J3>(Iyy3+Izz3))
	disp('WARNING: J3 must be <= Iyy3+Izz3')
      end
      disp(['Error in element properties number '... 
	    num2str(element(elnum).properties) ...
	    'used by element ' num2str(elnum) ' on line'...
	    num2str(element(elnum).lineno) '.'])
    end
    J=polyval(rn1*J1+rn2*J2+rn3*J3,rgpts(i));
    A=polyval(rn1*A1+rn2*A2+rn3*A3,rgpts(i));
    krod=krod+rgpw(i)*rodsfs*rodsfs'*A*E*Jac;
    ktor=ktor+rgpw(i)*rodsfs*rodsfs'*J*G*Jac;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % Derivation of Mass matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  numbeamgauss=numbeamgauss+3;
  [bgpts,bgpw]=gauss(numbeamgauss);
  mb1=zeros(6,6);
  % Local Bending in x-y plane
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1,bgpts(i));
             polyval(bn2,bgpts(i))*Jac;
             polyval(bn3,bgpts(i));
             polyval(bn4,bgpts(i))*Jac;
             polyval(bn5,bgpts(i));
             polyval(bn6,bgpts(i))*Jac];
    A=polyval(rn1*A1+rn2*A2+rn3*A3,bgpts(i));
    mb1=mb1+bgpw(i)*beamsfs*beamsfs'*rho*A*Jac;%pause
  end
  
  % Local Bending in x-z plane
  mb2=zeros(6,6);
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1,bgpts(i));
             -polyval(bn2,bgpts(i))*Jac;
             polyval(bn3,bgpts(i));
             -polyval(bn4,bgpts(i))*Jac;
             polyval(bn5,bgpts(i));
             -polyval(bn6,bgpts(i))*Jac];
    A=polyval(rn1*A1+rn2*A2+rn3*A3,bgpts(i));
    mb2=mb2+bgpw(i)*beamsfs*beamsfs'*rho*A*Jac;
  end
  
  % Local Extension in x, torsion about x
  numrodgauss=numrodgauss+1;
  [rgpts,rgpw]=gauss(numrodgauss);
  mrod=zeros(3,3);
  mtor=zeros(3,3);
  for i=1:numrodgauss
    rodsfs=[polyval(rn1,rgpts(i));
            polyval(rn2,rgpts(i));
            polyval(rn3,rgpts(i))];
    J=polyval(rn1*(Iyy1+Izz1)+rn2*(Iyy2+Izz2)+rn3*(Iyy3+Izz3),rgpts(i));
    A=polyval(rn1*A1+rn2*A2+rn3*A3,rgpts(i));
    mrod=mrod+rgpw(i)*rodsfs*rodsfs'*A*rho*Jac;
    mtor=mtor+rgpw(i)*rodsfs*rodsfs'*J*rho*Jac;
  end
  
  % Assembling each stiffness matrix into the complete elemental 
  % stiffness matrix
  k=zeros(18,18);
  k([2 6 8 12 14 18],[2 6 8 12 14 18])=kb1;
  k([3 5 9 11 15 17],[3 5 9 11 15 17])=kb2;
  k([1 7 13],[1 7 13])=krod;
  k([4 10 16],[4 10 16])=ktor;
  
  % Assembling each mass matrix into the complete elemental 
  % mass matrix
  m=zeros(18,18);
  m([2 6 8 12 14 18],[2 6 8 12 14 18])=mb1;
  m([3 5 9 11 15 17],[3 5 9 11 15 17])=mb2;
  m([1 7 13],[1 7 13])=mrod;
  m([4 10 16],[4 10 16])=mtor;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Coordinate rotations
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  R1=([x2 y2 z2]-[x1 y1 z1]);
  lam1=R1/norm(R1);
  R2=([x4 y4 z4]-[x1 y1 z1]);
  R2perp=R2-dot(R2,lam1)*lam1;
  udirec=0;
  while norm(R2perp)<10*eps
    udirec=udirec+1;
    %disp('oops')
    %pause
    [minval,minloc]=min(lam1);
    R2perp=zeros(1,3);
    R2perp(udirec)=1;
    R2perp=R2perp-dot(R2perp,lam1)*lam1;
  end
  lam2=R2perp/norm(R2perp);
  lam3=cross(lam1,lam2);
  lamloc=[lam1;lam2;lam3];
  lam=sparse(18,18);
  lam(1:3,1:3)=lamloc;
  lam(4:6,4:6)=lamloc;
  lam(7:9,7:9)=lamloc;
  lam(10:12,10:12)=lamloc;
  lam(13:15,13:15)=lamloc;
  lam(16:18,16:18)=lamloc;
  
% $$$     lam=[lamloc z z z z z;
% $$$          z lamloc z z z z;
% $$$          z z lamloc z z z;
% $$$          z z z lamloc z z;
% $$$          z z z z lamloc z;
% $$$          z z z z z lamloc];
  element(elnum).lambda=lam;
  element(elnum).m=m;
  element(elnum).k=k;

  kg=lam'*k*lam;
  mg=lam'*m*lam;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  bn1=bnodes(1);bn2=bnodes(2);bn3=bnodes(3);
  indices=[bn1*6+(-5:0) bn2*6+(-5:0) bn3*6+(-5:0)] ;


  K(indices,indices)=K(indices,indices)+kg;
  M(indices,indices)=M(indices,indices)+mg;

  % At this point we also know how to draw the element (what lines
  % and surfaces exist). For the beam3 element, 2 lines are
  % appropriate.
  numlines=size(lines,1);
  lines(numlines+1,:)=[bn1 bn3];
  lines(numlines+2,:)=[bn3 bn2];
  
%diag(M)
elseif strcmp(mode,'istrainforces')
  %disp('beam3 initial strain forces')
  % We need to have the stiffness matrix and the coordinate roation matrix.
  
  k=element(elnum).k;
  lam=element(elnum).lambda;
  kg=lam'*k*lam;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Add initial deflections
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if exist('distype')==0
    distype=0;
  end
  %disp('before loop')
  %bprops
  %length(bprops)
  if length(bprops)>15|length(bprops)==14
    %disp('in loop')
    if distype==0
      %disp('need to check for limits, then apply them/make them')
      %normal distribution
      if ~exist('ldx2')
        ldx2=inf;
        ldy2=inf;
        ldz2=inf;
        ldrx2=inf;
        ldry2=inf;
        ldrz2=inf;
        ldx3=inf;
        ldy3=inf;
        ldz3=inf;
        ldrx3=inf;
        ldry3=inf;
        ldrz3=inf;
      end
      dev=randn(1)*sdx2;
      while dev>ldx2
        dev=randn(1)*sdx2;
      end
      dx2=mdx2+dev;
      dev=randn(1)*sdy2;
      while dev>ldy2
        dev=randn(1)*sdy2;
      end
      dy2=mdy2+dev;
      dev=randn(1)*sdz2;
      while dev>ldz2
        dev=randn(1)*sdz2;
      end
      dz2=mdz2+dev;
      dev=randn(1)*sdrx2;
      while dev>ldrx2
        dev=randn(1)*sdrx2;
      end
      drx2=mdrx2+dev;
      dev=randn(1)*sdry2;
      while dev>ldry2
        dev=randn(1)*sdry2;
      end
      dry2=mdry2+dev;
      dev=randn(1)*sdrz2;
      while dev>ldrz2
        dev=randn(1)*sdrz2;
      end
      drz2=mdrz2+dev;
    elseif abs(distype)==1
      %uniform distribution
      %disp('uniform distribution')
      dx2=mdx2+(rand(1)-.5)*sdx2*2;%pause
      dy2=mdy2+(rand(1)-.5)*sdy2*2;
      dz2=mdz2+(rand(1)-.5)*sdz2*2;
      drx2=mdrx2+(rand(1)-.5)*sdrx2*2;
      dry2=mdry2+(rand(1)-.5)*sdry2*2;
      drz2=mdrz2+(rand(1)-.5)*sdrz2*2;
    else
      disp(['Invalid distribution type ' num2str(distype) '.'])
    end
    if exist('mdx3')==0
      %Calculate 'smooth' displacements of node 3 if they are
      %undefined

      node3disps=-k(13:18,13:18)\k(13:18,7:12)*[dx2;dy2;dz2;drx2;dry2;drz2];
      dx3=node3disps(1);
      dy3=node3disps(2);
      dz3=node3disps(3);
      drx3=node3disps(4);
      dry3=node3disps(5);
      drz3=node3disps(6);
    else
      if distype==0
        %normal distribution
        dev=randn(1)*sdx3;
        while dev>ldx3
          dev=randn(1)*sdx3;
        end
        dx3=mdx3+dev;
        dev=randn(1)*sdy3;
        while dev>ldy3
          dev=randn(1)*sdy3;
        end
        dy3=mdy3+dev;
        dev=randn(1)*sdz3;
        while dev>ldz3
          dev=randn(1)*sdz3;
        end
        dz3=mdz3+dev
        dev=randn(1)*sdrx3;
        while dev>ldrx3
          dev=randn(1)*sdrx3;
        end
        drx3=mdrx3+dev;
        dev=randn(1)*sdry3;
        while dev>ldry3
          dev=randn(1)*sdry3;
        end
        dry3=mdry3+dev;
        dev=randn(1)*sdrz3;
        while dev>ldrz3
          dev=randn(1)*sdrz3;
        end
        drz3=mdrz3+dev;
      elseif abs(distype)==1
        %uniform distribution
        dx3=mdx3+(rand(1)-.5)*sdx3*2;
        dy3=mdy3+(rand(1)-.5)*sdy3*2;
        dz3=mdz3+(rand(1)-.5)*sdz3*2;
        drx3=mdrx3+(rand(1)-.5)*sdrx3*2;
        dry3=mdry3+(rand(1)-.5)*sdry3*2;
        drz3=mdrz3+(rand(1)-.5)*sdrz3*2;
      else
        disp(['Invalid distribution type ' num2str(distype) '.'] )    
      end    
    end
    %disp('hello'),pause
    Fistrain=k*[0;0;0;0;0;0;dx2;dy2;dz2;drx2;dry2;drz2;dx3;dy3;dz3;drx3;dry3;drz3];
    % Need to transform Fistrain into global coordinates
    Fistrain=lam'*Fistrain;
    %  mg=lam'*m*lam;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Assembling matrices into global matrices
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bn1=bnodes(1);bn2=bnodes(2);bn3=bnodes(3);
    indices=[bn1*6+(-5:0) bn2*6+(-5:0) bn3*6+(-5:0)] ;
    
    if length(Fepsn)<max(indices);
      Fepsn(max(indices))=0;
    end
    Fepsn(indices)=Fepsn(indices)+Fistrain;
  end
  
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end
