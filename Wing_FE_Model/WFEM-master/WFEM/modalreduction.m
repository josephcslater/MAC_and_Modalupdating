function [Mred,Kred,T]=modalreduction(M,K,modes)
%[Mr,Kr,T]=MODALREDUCTION(M,K,modes) keeps only the modes listed in
%MODES. 
%[Mr,Kr,T]=MODALREDUCTION(M,K,nmodes) keeps only the first NMODES modes. 
% Mr=transpose(T)*M*T, Kr=transpose(T)*K*T  
% Reduced coordinate system forces can be obtained by 
% fr=transpose(T)*F
%
% Reduced damping matrix can be obtained using Cr=transpose(T)*C*T.
%
% If mode shapes are obtained for the reduced system, full system mode shapes 
% are phi=T*phi_r
%

%
% Copyright Joseph C. Slater, 6/12/2002
  nm=length(modes);
  if nm==1
    modes=1:modes;
  end
  
  ndof=length(M);
  % From $$$$$$ to $$$$$$ is the actual eigensolution.
  OPTS.issym=1;
  OPTS.isreal=1;
  %To increase accuracy: (Agnes)
  OPTS.tol=eps/10;
  [minvals,minvallocs]=sort(diag(K)./diag(M));
  
  shift=minvals(min(7,length(minvals)));
  
  [phi,d]=eigs((K+K')/2+shift*(M+M')/2,(M+M')/2,min([ size(K,1) ...
		    max([floor(sqrt(size(K,1))) 100])]),0,OPTS);
  d=d-shift;
  % $$$$$$
  %[phi,d]=eig(M\K);
  
  [S,index]=sort(diag(d));
  min(index)
  max(index)
  size(phi)
  
  size(modes)
  %index
  %size(phi)
  index
  modes
  disp('phi')
  phi=phi(:,index); 
  disp('T')
  T=phi(:,modes);
  size(T)
  Mred=T'*M*T;
  Kred=T'*K*T;

