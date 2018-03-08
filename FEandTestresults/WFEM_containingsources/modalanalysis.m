function [fs,fms]=modalanalysis
% MODALANALYSIS solves for the modes of the structure and plots
% them to the screen with a pretty lame menu. More later.
  
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
global F     % Prescribed forces
global Fnln  % Forces due to changing coordinate systems. 
global lines
global surfs
global bcs   % Boundary Conditions
global slavedofs
global fs
global fms 
global filename
global Tr 
global Mr
global Kr


i=1:size(K,1)
lslave=length(slavedofs)

if isempty(fs)|restart==0,% This is old stuff from previous restart
                          % code. If we haven't any natural
                          % frequencies already calculated, or we
                          % don't want a restart.
  if lslave>0% Here we need to reduce out constrained DOFs before
             % static analysis. Slave DOFs are defined when we apply
             % constraints or rigid elements.
    i(slavedofs)=zeros(length(slavedofs),1);%pop zeros into the k
                                            %indices corresponding
                                            %to slave dofs
    i=sort(i); % sort them. This puts the slave coordinates all up
               % top. 
%master1=[3;9;15;21;27;33;39;45;51;57;63;69;75;81;87;93;99;105;111;117;123;129;135;141;147;153;159;165;171;177;183;189;195;201;207;213;219;225;231;237;243;249]';
%length(master1)
%[1;2;3;7;8;9;13;14;15;19;20;21;25;26;27;31;32;33;37;38;39;43;44;45;49;50;51;55;56;57;61;62;63;67;68;69;73;74;75;79;80;81;85;86;87;91;92;93;97;98;99;103;104;105;109;110;111;115;116;117;121;122;123;127;128;129;133;134;135;139;140;141;145;146;147;151;152;153;157;158;159;163;164;165;169;171;171;175;176;177;181;182;183;187;188;189;193;194;195;199;200;201;205;206;207;211;212;213;217;218;219;223;224;225;229;230;231;235;236;237;241;242;243;247;248;249]'
%master2=[1:8];
master1=i(lslave+1:size(nodes,1)*6)%The master coordinates are
master2=[1:7];                                       %the ones that we didn't
                                       %reduce out via some
                                       %constraints. The master
                                       %coordinates are all of the
                                       %coordinates after the zero
                                       %values in i, but not those
                                       %beyond coordinates
                                       %corresponding to nodal
                                       %DOFS. The upper limit of
                                       %size(nodes,1)*6 fixes this
                                       %at the last real DOF.
   [Mr,Kr,T,master,slave]=guyan(M,K,master1)
[Mred,Kred,T,master,slave]=serep(Mr,Kr,master2)
  else
    [Mr,Kr,T,master,slave]=guyan(M,K)
    
  end
  % From $$$$$$ to $$$$$$ is the actual eigensolution.
  OPTS.issym=1;
  OPTS.isreal=1;
  %To increase accuracy: (Agnes)
  %OPTS.tol=eps/10;
%  [T,Mred,Kred]=reducedofs;
  [minvals,minvallocs]=sort(diag(Kred)./diag(Mred));
  
  shift=minvals(min(7,length(minvals)));
  [fms,f]=eigs((Kred+Kred')/2+shift*(Mred+Mred')/2,(Mred+Mred')/2,min([ size(Kred,1) ...
		    max([floor(sqrt(size(Kred,1))) 100])]),0,OPTS);
  f
  fs=sqrt(diag(f)-shift)/2/pi
  % $$$$$$
  fms=T*fms;
  [fs,fsindex]=sort(fs);
  fs
  fms=fms(:,fsindex);
else
  disp('Using previous eigensolution.')
end


