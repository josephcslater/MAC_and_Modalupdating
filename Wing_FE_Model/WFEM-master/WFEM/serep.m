function [Mred,Kred,T,master,slave]=serep(M,K,master)
%GUYAN Reduce size of second order system of equations by 
% applying Guyan reduction.
% M x'' + K x = 0
% Mr xm'' + Kr xm = 0
% Where x = T*xm, Mr= transpose(T)*M*T, Kr= transpose(T)*K*T
%
%[Mr,Kr,T,master,slave]=SEREP(M,K,master) slaves coordinates which are not
% listed in the vector master. Master coordinates are returned
% in the vector master, slave coordinate indices are returned in the vector
% slave.
%

%[Kr,T,master,slave]=SEREP(K,master) performs static condensation, slaving 
% coordinates that are not listed in the vector master. Master coordinates are 
% returned in the vector master, slave coordinate indices are returned in the 
% vector slave.
%
%[Mr,Kr,T,master,slave]=SEREP(M,K) slaves coordinates i for which 
% 0.1>(k(i,i)/m(i,i))/max(k(i,i)/m(i,i)). Master coordinates are returned
% in the vector master, slave coordinate indices are returned in the vector
% slave.
%
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
  nm=length(master);
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
  nm
  [S,index]=sort(diag(d));
  %index
  %size(phi)
  phi=phi(:,index);
  phitr=phi(nm+1:ndof,1:nm);
  phirr=phi(1:nm,1:nm);
  %disp('Eigensolution can be sped up')
  
  slave=(1:ndof)';
  slave(master)=zeros(nm,1);
  slave=sort(slave);
  slave=slave(nm+1:ndof);
  
  T=zeros(ndof,nm);
  T(master,1:nm)=eye(nm);
  T(slave,1:nm)=phitr/phirr;
  Mred=T'*M*T;
  Kred=T'*K*T;
end

