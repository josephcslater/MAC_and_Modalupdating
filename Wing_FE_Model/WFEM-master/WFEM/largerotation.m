function [rigidloc,flexloc]=largerotation(disp,modes)
%[rigidloc,flexloc]=largerotation(disp,modes)  
% Given large rotations (participations of modes 4-6), determine the
% actual position of the system. RIGIDLOC is the location of the
% system presuming it's rigid. FLEXLOC is the location including flexibility
  
% global restart
% global K
% global Ks
% global M
% global nodes
% global beamprops
% global matprops
% global points
% global elprops
% global element
% global Fepsn % Forces representing initial strain
% global F     % Prescribed forces
% global Fnln  % Forces due to changing coordinate systems. 
% global lines
% global surfs
% global bcs   % Boundary Conditions
% global slavedofs
% global fs
% global fms 
% global filename
% global Tr 
% global Mr
% global Kr

thetax=modes(4,4)*modes(:,4)\disp;
thetay=modes(5,5)*modes(:,5)\disp;
thetaz=modes(6,6)*modes(:,6)\disp;

%thetax transformation matrix
txm=[1       0       0;
     0 cos(thetax) -sin(thetax);
     0 sin(thetax) cos(thetax)];

%thetay transformation matrix
tym=[cos(thetay)  0  sin(thetay);
     0            1       0;
     -sin(thetay) 0  cos(thetay)];

%thetaz transformation matrix
tyz=[cos(thetaz) -sin(thetaz)  0;
     sin(thetaz) cos(thetaz)   0;
     0 0 1];


tm=txm*tym*tzm; % This doesn't work if more than one angle is
                % non-zero
		
% Next task is to move those nodal starting points (include that
% first mode) then add the displacements (in the correced
% directions). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate initial 'shape' vector - rotated shape without
%deformations

ishape=emptymode;
ishape(1:6:snodes*6)=nodes(:,1);
ishape(2:6:snodes*6)=nodes(:,2);
ishape(3:6:snodes*6)=nodes(:,3);

% Now rotate it. 
Trotate=spalloc(snodes*6,snodes*6,36*snodes);
for i=1:snodes;
  Trotate(i*6+[-5:0],i*6+[-5:0])=tm;
end

rigidloc=Trotate*ishape;


%Now we have the nodal locations after rotating. We need the nodal
%locations after rotating PLUS deforming


flexloc=rigidloc+Trotate*disp;