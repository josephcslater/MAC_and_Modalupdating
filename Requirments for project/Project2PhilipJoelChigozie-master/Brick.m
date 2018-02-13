function out=Brick(mode,b,c,d,e)
  
% Brick does as listed below.
% Brick properties (bprops) are in the order
% bprops=[E G rho]

%%
% Defining beam element properties in wfem input file:
% element properties
%   E G rho 
%
% Defining Brick element in wfem input file:
%   node1 node2 node3 node4 node5 node6 node7 node8 materialnumber 

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
global surfs
%
% Variables (local):
% ------------------
% bnodes  :    node/point numbers for actual brick nodes 1-8
% k       :    stiffness matrix in local coordiates
% kg      :    stiffness matrix rotated into global coordinates
% m       :    mass matrix in local coordiates
% mg      :    mass matrix rotated into global coordinates
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out=0;
if strcmp(mode,'numofnodes')
    % This allows a code to find out how many nodes this element has
    out=8;
end
if strcmp(mode,'generate')
  elnum=c;%When this mode is called, the element number is the 3rd
          %argument. 
  
          %The second argument (b) is the element
          %definition. For this element b is
          %node1 node2 node3 node4 node5 node6 node7 node8 materialnumber
  
          %There have to be 9 elements for this element's
          %definition (above)
  if length(b)==9
      element(elnum).nodes=b(1:8);
      element(elnum).properties=b(9);
  else 
	  b
      %There have to be nine numbers on a line defining the
      %element. 
      warndlg(['Element ' num2str(elnum) ' on line ' ...
               num2str(element(elnum).lineno) ' entered incorrectly.'], ...
              ['Malformed Element'],'modal')
      return
  end
 
end

% Here we figure out what the beam properties mean. If you need
% them in a mode, that mode should be in the if on the next line.
if strcmp(mode,'make')||strcmp(mode,'istrainforces')
  elnum=b;% When this mode is called, the element number is given
          % as the second input.
bnodes=[element(elnum).nodes];% The point is
                                                     % referred to
                                                     % as node 4
                                                     % below,
                                                     % although it
                                                     % actually
                                                     % calls the
                                                     % array points
                                                     % to get its
                                                     % location. Its
                                                     % not really a
                                                     % node, but
                                                     % just a point
                                                     % that helps
                                                     % define
                                                     % orientation. Your
                                                     % element may
                                                     % not need
                                                     % such a
                                                     % reference point.
  bprops=elprops(element(elnum).properties).a;% element(elnum).properties 
                                              % stores the
                                              % properties number
                                              % of the current
                                              % elnum. elprops
                                              % contains this
                                              % data. This is
                                              % precisely the
                                              % material properties
                                              % line in an
                                              % array. You can pull
                                              % out any value you
                                              % need for your use. 
  
% 
  if length(bprops)==3
      E=bprops(1);
      G=bprops(2);
      rho=bprops(3);
  else
      warndlg(['The number of material properties set for ' ...
               'this element (' num2str(length(bprops)) ') isn''t ' ...
               'appropriate for a beam3 element. '      ...
               'Please refer to the manual.'],...
              'Bad element property definition.','modal');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brick properties (bprops) are in the order
% bprops=[E G rho]


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
  x4=nodes(bnodes(4),1);
  y4=nodes(bnodes(4),2);
  z4=nodes(bnodes(4),3);
  x5=nodes(bnodes(5),1);
  y5=nodes(bnodes(5),2);
  z5=nodes(bnodes(5),3);
  x6=nodes(bnodes(6),1);
  y6=nodes(bnodes(6),2);
  z6=nodes(bnodes(6),3);
  x7=nodes(bnodes(7),1);
  y7=nodes(bnodes(7),2);
  z7=nodes(bnodes(7),3);
  x8=nodes(bnodes(8),1);
  y8=nodes(bnodes(8),2);
  z8=nodes(bnodes(8),3);

  xvec=[x1 x2 x3 x4 x5 x6 x7 x8]; % stores nodal x values in J_brick format
  yvec=[y1 y2 y3 y4 y5 y6 y7 y8]; % stores nodal y values in J_brick format
  zvec=[z1 z2 z3 z4 z5 z6 z7 z8]; % stores nodal z values in J_brick format
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Define gauss points and loop to integrate for Ke and Me
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Emat=E_matrix(E,G); % Calculates the E_matrix for the element material
  Me=zeros(24,24); % Preallocates the element mass matrix
  Ke=zeros(33,33); % Preallocates the element stiffness matrix
  
  numbrickgauss=5; % Chooses a number of gauss points for the stiffness matrix
  [bgpts,bgpw]=gauss([numbrickgauss,numbrickgauss,numbrickgauss]); % finds the gauss points and weights
  
  J0=J_brick(0,0,0,xvec,yvec,zvec);
  for i=1:size(bgpts,1) % loops over all the gauss points
      [J,dNdx,dNdy,dNdz]=J_brick(bgpts(i,1),bgpts(i,2),bgpts(i,3),xvec,yvec,zvec); % Finds J for the current Gauss point
      B=B_brick(dNdx,dNdy,dNdz); % Finds B for the current Gauss point
      Bt=B'; % Calculates B transpose
      Ki=bgpw(i)*Bt*Emat*B*det(J); % Calculates the weighted Gauss point stiffness
      Ke(1:24,1:24)=Ke(1:24,1:24)+Ki(1:24,1:24); % adds the weighted Gauss point stiffness to the element stiffness

      Ba=Ba_brick(J,bgpts(i,1),bgpts(i,2),bgpts(i,3)); 
      B=[B Ba];
      Bt=B'; % Calculates B transpose
      Ki=bgpw(i)*Bt*Emat*B*det(J0); % Calculates the weighted Gauss point stiffness
      Ke(1:33,25:33)=Ke(1:33,25:33)+Ki(1:33,25:33); % adds the weighted Gauss point stiffness to the element stiffness
      Ke(25:33,1:24)=Ke(25:33,1:24)+Ki(25:33,1:24);
  end
  
  numbrickgauss=numbrickgauss+1; % Adds more gauss points for the mass matrix
  [bgpts,bgpw]=gauss([numbrickgauss,numbrickgauss,numbrickgauss]); % finds the gauss points and weights
  for i=1:size(bgpts,1) % loops over all the gauss points
      J=J_brick(bgpts(i,1),bgpts(i,2),bgpts(i,3),xvec,yvec,zvec); % Finds J for the current Gauss point
      N=N_brick(bgpts(i,1),bgpts(i,2),bgpts(i,3)); % Finds N for the current Gauss point
      Nt=N'; % Calculates N transpose
      Mi=bgpw(i)*Nt*rho*N*det(J); % Calculates the weighted Gauss point mass
      Me=Me+Mi; % adds the weighted Gauss point mass to the element mass
  end
  
  [Ke]=ZeroDOFs(Ke);
  [Ke,Me]=Guyan_Brick(Ke,Me);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  bn1=bnodes(1);
  bn2=bnodes(2);
  bn3=bnodes(3);
  bn4=bnodes(4);
  bn5=bnodes(5);
  bn6=bnodes(6);
  bn7=bnodes(7);
  bn8=bnodes(8);
  indices=[bn1*6+(-5:-3) bn2*6+(-5:-3) bn3*6+(-5:-3) bn4*6+(-5:-3)...
           bn5*6+(-5:-3) bn6*6+(-5:-3) bn7*6+(-5:-3) bn8*6+(-5:-3)] ;


  K(indices,indices)=K(indices,indices)+Ke;
  M(indices,indices)=M(indices,indices)+Me;

  % At this point we also know how to draw the element (what lines
  % and surfaces exist). For the beam3 element, 2 lines are
  % appropriate. Just add the pair of node numbers to the lines
  % array and that line will always be drawn.
  numlines=size(lines,1);
  lines(numlines+1,:)=[bn1 bn2];
  lines(numlines+2,:)=[bn2 bn3];
  lines(numlines+3,:)=[bn3 bn4];
  lines(numlines+4,:)=[bn4 bn1];
  lines(numlines+5,:)=[bn5 bn6];
  lines(numlines+6,:)=[bn6 bn7];
  lines(numlines+7,:)=[bn7 bn8];
  lines(numlines+8,:)=[bn8 bn5];
  lines(numlines+9,:)=[bn1 bn5];
  lines(numlines+10,:)=[bn2 bn6];
  lines(numlines+11,:)=[bn3 bn7];
  lines(numlines+12,:)=[bn4 bn8];
  
  %If I have 4 nodes that I want to use to represent a surface, I
  %do the following.
  panelcolor=[1 0 1];% This picks a color. You can change the
                     % numbes between 0 and 1. 
  %Don't like this color? Use colorui to pick another one. Another
  %option is that if we can't see the elements separately we can
  %chunk up x*y*z, divide by x*y*x of element, see if we get
  %integer powers or not to define colors that vary by panel. 
  
  
  % You need to uncomment this line and assign values to node1,
  % node2, node3, and node4 in order to draw A SINGLE SURFACE. For
  % a brick, you need 6 lines like this. 
  surfs=[surfs;bn1 bn2 bn3 bn4 panelcolor];
  surfs=[surfs;bn2 bn6 bn7 bn3 panelcolor];
  surfs=[surfs;bn6 bn5 bn8 bn7 panelcolor];
  surfs=[surfs;bn5 bn1 bn4 bn8 panelcolor];
  surfs=[surfs;bn4 bn8 bn7 bn3 panelcolor];
  surfs=[surfs;bn1 bn5 bn6 bn2 panelcolor];
  
  %Each surface can have a different color if you like. Just change
  %the last three numbers on the row corresponding to that
  %surface. 

%diag(M)
elseif strcmp(mode,'istrainforces')
  % You don't need this
  % We need to have the stiffness matrix and the coordinate roation matrix.
 

  
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end
