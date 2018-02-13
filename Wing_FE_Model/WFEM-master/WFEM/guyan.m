function [Mr,Kr,T,master,slave,zerodofs]=guyan(m,k,thresh)
%GUYAN Reduce size of second order system of equations by 
% applying Guyan reduction.
% M x'' + K x = 0
% Mr xm'' + Kr xm = 0
% Where x = T*xm, Mr= transpose(T)*M*T, Kr= transpose(T)*K*T
%
%[Mr,Kr,T,master,slave,zerodofs]=GUYAN(M,K,frac) slaves the
% fraction of coordinates above a to-be-determined
% thresh>(k(i,i)/m(i,i))/max(k(i,i)/m(i,i)).  Master coordinates are returned
% in the vector master, slave coordinate indices are returned in the vector
% slave.
%
%[Mr,Kr,T,master,slave,zerodofs]=GUYAN(M,K,master) slaves coordinates which are not
% listed in the vector master. Master coordinates are returned
% in the vector master, slave coordinate indices are returned in the vector
% slave.
%
%[Kr,T,master,slave,zerodofs]=GUYAN(K,master) performs static condensation, slaving 
% coordinates that are not listed in the vector master. Master coordinates are 
% returned in the vector master, slave coordinate indices are returned in the 
% vector slave.
%
%[Mr,Kr,T,master,slave,zerodofs]=GUYAN(M,K) removes only degrees of
% freedom for which the stiffnesses are zero. See below. This
% adustment is always taken.
%
% Reduced coordinate system forces can be obtained by 
% fr=transpose(T)*F
%
% Reduced damping matrix can be obtained using Cr=transpose(T)*C*T.
%
% If mode shapes are obtained for the reduced system, full system mode shapes 
% are phi=T*phi_r
%
% The output argument zerodofs is an option list of degrees of
% freedom that have zero stiffness. This happens when an element
% with 3DOF nodes is used. The other 3DOFS have no stiffness
% (eg. rotational coordinates in a rod). guyan will also strip out
% those DOFs from the model by eliminating the columns
% corresponding to those DOFS, and 

% Modified 9/8/2003 to ditch unconstrained coordinates (e.g. rod
% element issues)
% Copyright Joseph C. Slater, 6/19/2000

% Check if the second incoming argument is actually just a list of
% coordinates.   

  
if size(k,1)==1 | size(k,2)==1
  statcond=1;
  master=k;
  k=m;
  m=speye(size(k));
else
  statcond=0;
end


sprse=1;
%Make the matrices sparse, if they're not. Remember to expand them
%later if they came in full.
if ~issparse(m)
  m=sparse(m);
  k=sparse(k);
  sprse=0;
end

% If no reduction is ordered, remove only zero stiffness DOFs
if nargin==2
  thresh=0;
end

ncoord=length(m);


dm=diag(m);
dk=diag(k);
rat=dm./(dk+realmin);
mr=rat./max(rat);
mth=min(rat)/max(rat);


thresh(1);
autothresh=0;
if length(thresh)==1&thresh(1)~=0% Determine the master and slave
                                 % coordinates if som percentage,
                                 % thresh, are to be slaved. If
                                 % none are to be discarded, we
                                 % don't need to do much.
  if thresh>=1
    thresh=.9;%default slave rate
    warndlg({'thresh must be less than 1.';...
	     'thresh has been set to .9'; ...
	     'slaving 90 percent.'},'Threshold too high')
    %        disp('move on')
  autothresh=1;
  end

%  master=(mr)>thresh;
  numofmasters=floor((1-thresh)*ncoord);
  [rat,i]=sort(mr);
  slave=i(1:ncoord-numofmasters)';
  master=i(ncoord-numofmasters+1:ncoord)';
  lmaster=length(master);
elseif thresh(1)==0
  master=1:ncoord;lmaster=ncoord;
  slave=[];
else
  master=thresh;
  i=1:ncoord;
  lmaster=length(master);
  i(master)=zeros(lmaster,1);
  i=sort(i);
  slave=i(lmaster+1:ncoord);
end
% Find zero stiffness DOFs
%zerovals=full((diag(k)==0))';
zerovals=full(sum(abs(k))==0);
%zerovals
numzeros=sum(zerovals); % Number of on-stiffened DOFs
i=1:ncoord;

zerolocs=sort(i.*zerovals);
zerolocs=zerolocs((ncoord+1-numzeros):ncoord);
%zerolocs=[];
% checking should be added such that a constrained DOF dos not also
% get eliminated.Perhaps intersect with slave and don't mess with
% slaved coordinates. 
%disp(['line 131 of guyan, no checking of intersection of 0 stiffness' ...
%      ' DOFS and comparing against constrained DOFS.'])

if autothresh==1
  figure(gcf)
  plot((1:ncoord)',rat./max(rat),[1 ncoord],[thresh thresh])
  axis([1,ncoord, 0,1])
  legend('Normalized Ratios','Cutoff for Slave Coordinates',0)
  grid on
  xlabel('Coordinate number')
  ylabel('Normalized diagonal ratio.')
%   h=warndlg({['thresh must be greater than ' num2str(mth) ' for']...
%              ;'this problem. Run aborted.'},'Threshold too hign');
%   t=(1:.1:70)';
%   noise=sin(10*t)+sin(3*t);
%   quiet=0*t;
%   sound([noise;quiet;noise;quiet;noise])
%  slave=0;master=0;Mr=0;Kr=0;T=0;
%  pause
  %break
end
%T=[];
%T=sparse(T);
%T=zeros(ncoord,lmaster);
%master
master2=union(zerolocs,master);
slave2=setdiff(slave,zerolocs);
master=master2;
slave=slave2;
lmaster=length(master);
kmm=k(master,master);
ksm=k(slave,master);
kss=k(slave,slave);

kms=ksm';
%lmaster
%ncoord
if lmaster==ncoord
  T=speye(ncoord);
else
  T=sparse(length(master)+length(slave),length(master));%spy(T);pause
  T(master,1:lmaster)=speye(lmaster,lmaster);%spy(T);pause
  T(slave,1:lmaster)=-kss\ksm;%spy(T);pause
end


% Need to eliminate columns and zero the corresponding rows for
% null rotational DOFS.
% Zeroing rows.

size(sparse(length(zerolocs),length(master)));
T(zerolocs,1:lmaster)=sparse(length(zerolocs),length(master));
% removing columns
%Need to move the differences from the slave cordinates to the
%master coordinates.
if length(intersect(full(master),full(zerolocs)))<length(full(zerolocs))
  t=(1:.1:70)';
  noise=sin(10*t)+sin(3*t);
  quiet=0*t;
  sound([noise;quiet;noise;quiet;noise])
  disp('Warning, a constrained DOF is also zero-stiffness.')
  pause
end
ismasterzeros=ismember(master,zerolocs);
i=1:length(master);
masterkeepers=sort(abs(ismasterzeros-1).*i);

masterkeepers=masterkeepers(sum(ismasterzeros)+1:length(master));
T=T(:,masterkeepers);

if statcond~=1
  Mr=T'*m*T;
end
Kr=T'*k*T;
if sprse==0
  Mr=full(Mr);
  Kr=full(Kr);
  T=full(T);
end
if statcond==1
  Mr=Kr;
  Kr=T;
  T=master;
  master=slave;
end
master=setdiff(master,zerolocs);
zerodofs=zerolocs;
