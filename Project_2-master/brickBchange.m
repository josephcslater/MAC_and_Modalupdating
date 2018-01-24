function out=brickBchange(mode,b,c,d,e)
  
% brick does as listed below.  
% brick properties (bprops) are in the order
% bprops=[E G rho ]
%%
% Defining beam element properties in wfem input file:
% element properties
%   E G rho 
%
% Defining brick element in wfem input file:
%   node1 node2 node3 node4 node5 node6 node7 node8 and material#
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
% bnodes  :    node/point numbers for actual beam nodes 1-8 and material#
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
    % This allows a code to find out how many nodes this element has
    out=8;
end
if strcmp(mode,'generate')
  elnum=c;%When this mode is called, the element number is the 3rd
          %argument. 
  
          %The second argument (b) is the element
          %definition. For this element b is
          %node1 node2 node3 node4 node5 node6 node7 node8 and material#
  
          %There have to be 9 numbers for this element's
          %definition (above)
    if length(b)==9
      element(elnum).nodes=b(1:8); 
      element(elnum).properties=b(9);
    end
 
end

% Here we figure out what the brick properties mean. If you need
% them in a mode, that mode should be in the if on the next line.
if strcmp(mode,'make')||strcmp(mode,'istrainforces')
  elnum=b;% When this mode is called, the element number is given
          % as the second input.
  bnodes=[element(elnum).nodes];
  bprops=elprops(element(elnum).properties).a;

E=bprops(1);
G=bprops(2);
rho=bprops(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode,'make')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Define beam node locations for easy later referencing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   numofnodes = 8;
  x = [];
  y = [];
  z = [];
  for i=1:numofnodes   
   xi = nodes(bnodes(i),1);
   yi = nodes(bnodes(i),2);
   zi = nodes(bnodes(i),3);
   x = [x xi];%x locations in cartesian coords
   y = [y yi];%y locations in cartesian coords
   z = [z zi];%z locations in cartesian coords
  end
  x = x';
  y = y';
  z = z';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %E matrix formulation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  v = (E/(2*G))-1; %Poisson's Ratio
  v1 = 1-v;
  v2 = 1+v;
  v3 = 1-2*v;
  
  Ei = zeros(6,6);
  Ei(1:3,1:3) = v;
  i = 1;
  for i = 1:3
  Ei(i,i) = v1;
  Ei((i+3),(i+3)) = v3/2;
  end
  Ematrix = Ei*E/(v2*v3);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Determine Mass & Stiffness Matricies using Gauss Integration
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  gaussnum = 5; % Number of Gauss points for integration of brick element
  [bgpts,bgpw] = gauss([gaussnum,gaussnum,gaussnum]);
  Me=zeros(24,24);
  Ke=zeros(24,24);
 
  Xi =   [-1 -1 -1 -1  1  1  1 1];
  Etai = [-1 -1  1  1 -1 -1  1 1];
  Zi =   [ 1 -1 -1  1  1 -1 -1 1];
  
  for i=1:size(bgpts,1)  %loop for gauss integration
      
      Xe = bgpts(i,1);
      Eta = bgpts(i,2);
      Z = bgpts(i,3);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Determine N Matrix
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      N = []; 
      n = 1;
      for n = 1:numofnodes
          Ni = 0.125*(1+Xe*Xi(n))*(1+Eta*Etai(n))*(1+Z*Zi(n)); 
          Ni = [Ni Ni Ni];
          N = [N diag(Ni)];
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Determine J Matrix
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      dNdX = [];
      dNdN = [];
      dNdZ = [];
      n = 1;
      for n = 1:numofnodes
          dNdXi = 0.125*Xi(n)*(1+Eta*Etai(n))*(1+Z*Zi(n));
          dNdNi = 0.125*Etai(n)*(1+Xe*Xi(n))*(1+Z*Zi(n));
          dNdZi = 0.125*Zi(n)*(1+Xe*Xi(n))*(1+Eta*Etai(n));
          dNdX = [dNdX dNdXi];
          dNdN = [dNdN dNdNi];
          dNdZ = [dNdZ dNdZi];
      end
      J = [dNdX; dNdN; dNdZ]*[x y z];
      J_inverse = J^-1;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Determine B Matrix
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Alpha=zeros(6,9);
     Alpha(1,1)=1;Alpha(2,4)=1;Alpha(3,7)=1;Alpha(4,2)=1;Alpha(4,4)=1;Alpha(5,6)=1;Alpha(5,8)=1;Alpha(6,3)=1;Alpha(6,7)=1;
     zero=zeros(3,3);
     Beta=[J_inverse zero zero;zero J_inverse zero;zero zero J_inverse];
     Gamma=zeros(9,24);
     Gamma(1,1)=dNdX(1);Gamma(4,2)=dNdX(1);Gamma(7,3)=dNdX(1);
     Gamma(1,4)=dNdX(2);Gamma(4,5)=dNdX(2);Gamma(7,6)=dNdX(2);
     Gamma(1,7)=dNdX(3);Gamma(4,8)=dNdX(3);Gamma(7,9)=dNdX(3);
     Gamma(1,10)=dNdX(4);Gamma(4,11)=dNdX(4);Gamma(7,12)=dNdX(4);
     Gamma(1,13)=dNdX(5);Gamma(4,14)=dNdX(5);Gamma(7,15)=dNdX(5);
     Gamma(1,16)=dNdX(6);Gamma(4,17)=dNdX(6);Gamma(7,18)=dNdX(6);
     Gamma(1,19)=dNdX(7);Gamma(4,20)=dNdX(7);Gamma(7,21)=dNdX(7);
     Gamma(1,22)=dNdX(8);Gamma(4,23)=dNdX(8);Gamma(7,24)=dNdX(8);
     Gamma(2,1)=dNdN(1);Gamma(5,2)=dNdN(1);Gamma(8,3)=dNdN(1);
     Gamma(2,4)=dNdN(2);Gamma(5,5)=dNdN(2);Gamma(8,6)=dNdN(2);
     Gamma(2,7)=dNdN(3);Gamma(5,8)=dNdN(3);Gamma(8,9)=dNdN(3);
     Gamma(2,10)=dNdN(4);Gamma(5,11)=dNdN(4);Gamma(8,12)=dNdN(4);
     Gamma(2,13)=dNdN(5);Gamma(5,14)=dNdN(5);Gamma(8,15)=dNdN(5);
     Gamma(2,16)=dNdN(6);Gamma(5,17)=dNdN(6);Gamma(8,18)=dNdN(6);
     Gamma(2,19)=dNdN(7);Gamma(5,20)=dNdN(7);Gamma(8,21)=dNdN(7);
     Gamma(2,22)=dNdN(8);Gamma(5,23)=dNdN(8);Gamma(8,24)=dNdN(8);
     Gamma(3,1)=dNdZ(1);Gamma(6,2)=dNdZ(1);Gamma(9,3)=dNdZ(1);
     Gamma(3,4)=dNdZ(2);Gamma(6,5)=dNdZ(2);Gamma(9,6)=dNdZ(2);
     Gamma(3,7)=dNdZ(3);Gamma(6,8)=dNdZ(3);Gamma(9,9)=dNdZ(3);
     Gamma(3,10)=dNdZ(4);Gamma(6,11)=dNdZ(4);Gamma(9,12)=dNdZ(4);
     Gamma(3,13)=dNdZ(5);Gamma(6,14)=dNdZ(5);Gamma(9,15)=dNdZ(5);
     Gamma(3,16)=dNdZ(6);Gamma(6,17)=dNdZ(6);Gamma(9,18)=dNdZ(6);
     Gamma(3,19)=dNdZ(7);Gamma(6,20)=dNdZ(7);Gamma(9,21)=dNdZ(7);
     Gamma(3,22)=dNdZ(8);Gamma(6,23)=dNdZ(8);Gamma(9,24)=dNdZ(8);
     B=Alpha*Beta*Gamma;
     
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Stiffness
      Bt=B';
      Ki = bgpw(i)*Bt*Ematrix*B*det(J);
      Ke = Ke+Ki;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Mass 
      Nt=N';
      Mi = bgpw(i)*Nt*rho*N*det(J);
      Me = Me+Mi;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Assembling matrices into global matrices
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


% Me
% Ke
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
  %surfs=[surfs;node1 node2 node3 node4 panelcolor];
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
