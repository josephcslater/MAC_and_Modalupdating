function [thetas,fildisp]=filterrbmodes(disp,modes)
%[thetas,disp]=filterrbmodes(disp,modes)
% Given large rotations (participations of modes 4-6), determine the
% deflection of the system relative to the rigid body
% location. thetas is the location of the
% system presuming its rigid. disp is th displacement vector
  
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
theta=[thetax;thetay;thetaz];
  
fildisp=disp-modes(:,4)*(modes(:,4)\disp);
