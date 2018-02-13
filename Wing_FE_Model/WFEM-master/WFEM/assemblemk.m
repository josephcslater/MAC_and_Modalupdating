function assemblemk
% This simply assembles the mass and stiffness matrices. 
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
global Cd Cv Ca nout Bso nin
global Bso
global Bthermal
global nout

filename;

if exist('restartsamemk')==0
  restartsamemk=0;
end
%size(K),pause

%filenameroot=filename(1:findstr(filename,'.')-1),%pause
%if exist([filenameroot '.mat'])==2
%  sourcefileinfo=dir(filename);
%  savefileinfo=dir([filenameroot '.mat']);
%else
%  savefileinfo=dir(filename);
%  sourcefileinfo=savefileinfo;
%end

%if size(K,1)==0|restartsamemk==0|datenum(savefileinfo.date)< ...
%      datenum(sourcefileinfo.date)
%  if datenum(savefileinfo.date)<datenum(sourcefileinfo.date)
%    disp(['Warning: M and K may have changed, but I''ll reuse them ' ...
%          'anyway'])
%    beep, beep, beep
%  end
%  


  K=sparse(size(nodes,1)*6,size(nodes,1)*6);
  M=K;
  Ca=sparse(0,size(nodes,1)*6);
  Cd=Ca;
  Cv=Ca;
  nout=0;
  Bso=Ca';
  Bthermal=Bso;
  nout=0;
  nin=0;
  numelements=length(element);
  waitbar(0,'Assembling Elements')
  for elnum=1:numelements
    h=waitbar(elnum/numelements);
    fcall=[element(elnum).elementtype '(''make'',' num2str(elnum) ');'];
%     if findstr(fcall,'panel')>0
%       fcall,pause
%     end
    
    eval(fcall);
  end
  delete(h)
  %  eval(['save ' filename(1:findstr(filename,'.')-1) '.mat'])
%else
  %eval(['load ' filename(1:findstr(filename,'.')-1) '.mat M K;'])
%  disp('Reading previous M and K.')%,pause
%end
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
