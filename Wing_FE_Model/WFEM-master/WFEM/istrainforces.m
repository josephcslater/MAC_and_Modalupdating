function istrainforces
% This generates the forces needed to have initial strains.
global restart
global K
global Ks
global M
global nodes
global beamprops
global matprops
global points
global elprops
global element
global Fepsn % Forces representing initial strain
global F     % Prescribe forces
global Fnln  % Forces due to changing coordinate systems. 
global lines
global surfs
global bcs %Boundary Conditions 
global filename
global restartsamemk

Fepsn=sparse(size(K,1),1);
%if length(Fepsn)<size(K,1)
%  Fepsn(size(K,1),1)=0;
%end
%disp('hello'),pause
numelements=length(element);
%waitbar(0,'Generating Initial Strain Forces')
for elnum=1:numelements
 % h=waitbar(elnum/numelements);
    fcall=[element(elnum).elementtype '(''istrainforces'',' ...
           num2str(elnum) ');'];
    eval(fcall);
    %pause
end

%delete(h)
sk=size(K,1);
if length(F)<sk
  F(sk)=0;
end
if length(Fepsn)<sk
  Fepsn(sk)=0;
end
if length(Fnln)<sk
  Fnln(sk)=0;
end
