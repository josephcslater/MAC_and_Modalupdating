function applyconstraint(a,b,c)
% APPLYCONSTRAINT applies constraint conditions between two nodes using the Lagrange
% Multiplier technique, see Cook, Malkus and Plesha. See wfem.pdf for details. 

% Copyright Joseph C. Slater, 2002
  
% We're going to act on a number of these, so we need access to
% them. K will be appended to with the boundary condition
% augmentation terms. The other matrices need to have zeros
% appended to keep them of the correct dimension. It's imperative
% to understand that static condensation MUST be applied for
% dynamic analysis to be performed. 
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
global ccs   % Constraint Conditions
global slavedofs
global fido
global Bthermal

if strcmp('read',a)
  % This is going to read in all of the constraint condition data.
  command=b;
  fid=c;
  line=fgetl(fid);
  numberofnumbers=length(str2num(line));
  i=0;
  % We're loading all of the data in at once. We'll end up with
  % chunks of like constraint conditions, separated only if the user
  % does so. I'd like to separate all of these later if
  % possible. This will allow spiffier plotting.
  while strcmp(line(1),'%')
    % We can't increment i here because it is used for placing the
    % data (see below)
    line=fgetl(fid);
  end
  numline=eval(['[' line '];' ]);
  while length(numline)==numberofnumbers&(numline~=-1)% the length
                                                      % assures
                                                      % that we're
                                                      % continuing
                                                      % to get
                                                      % similar
                                                      % type of date
    i=i+1;
    ccdata(i,:)=str2num(line);
    line=fgetl(fid);
    if ~isstr(line)%EOL returns a numerical -1
      line=num2str(line);
    end
    while length(line)>0&&strcmp(line(1),'%')
      line=fgetl(fid);
      if ~isstr(line)
        line=num2str(line);
      end
    end
    numline=eval([ '[' line '];' ]);
  end

  if exist('ccs')
    ccnum=length(ccs)+1;
  else
    ccnum=1;
  end
  
  if length(findstr(b,'clamp'))==1
    ccs(ccnum).type='clamp';
    ccs(ccnum).ccdata=ccdata;
  elseif length(findstr(b,'pin'))==1
    ccs(ccnum).type='pin';
    ccs(ccnum).ccdata=ccdata;
  elseif length(findstr(b,'roller'))==1
    ccs(ccnum).type='roller';
    ccs(ccnum).ccdata=ccdata;
  elseif length(findstr(b,'ball'))==1
    ccs(ccnum).type='ball';
    ccs(ccnum).ccdata=ccdata;
  elseif length(findstr(b,'surfaceball'))==1
    ccs(ccnum).type='surfaceball';
    ccs(ccnum).ccdata=ccdata;
  elseif length(findstr(b,'surface'))==1
    ccs(ccnum).type='surface';
    ccs(ccnum).ccdata=ccdata;
  elseif length(findstr(b,'rbeam'))==1
    ccs(ccnum).type='rbeam';
    ccs(ccnum).ccdata=ccdata;
  elseif length(findstr(b,'rod'))==1
    ccs(ccnum).type='rod';
    ccs(ccnum).ccdata=ccdata;
  elseif length(findstr(b,'rigidbody'))==1
    ccs(1).type='rigidbody'; % Since we won't use any other
                             % conditions, why care about them at
                             % all. 
    ccs(1).ccdata=ccdata;
  end
elseif strcmp('apply',a)
  if exist('slavedofs')==0
    slavedofs=[];
  end
  % This is where we apply all of the boundary conditions. 
  for i=1:length(ccs)
    ccnumb=i;% This is the index number of the
             % boundary condition under consideration
    b=ccs(i).type;
    data=ccs(i).ccdata;
    numccs=size(data,1); %Number of boundary conditions
    if length(findstr(b,'rigidbody'))==1
      %If we have a rigid body, no other constraint need to be
      %applied. What we do here is apply a rigid beam constrain
      %between each and every node and the 'primary' node as
      %defined by data(1). Note the similarity of this loop to that
      %for a rigid beam connection. 
      centernode=data(1,1);

      numccs=size(nodes,1);
      %need to exclude the node we are fixing to from the list of ...
	%  nodes. 
      listofnodes=1:size(nodes,1);
      
      discardednodes=sort((listofnodes~=centernode).*listofnodes);%,pause
      
      	for j=discardednodes(2:numccs)
	  sk=size(K,1);
	  nodenum1=centernode;
	  nodenum2=j;
	  tudv=nodes(nodenum2,:)-nodes(nodenum1,:);
	  % tangent unit direction vector is a bad variable name, but it's not being used elsewhere ...
	  rx=tudv(1);
	  ry=tudv(2);
	  rz=tudv(3);
	  skewsym=[0   rz  -ry;
		   -rz 0    rx;
		   ry  -rx   0];
	  conmat=[eye(6,6) [-eye(3,3) skewsym;zeros(3,3) ...
			    -eye(3,3)]];
	  K([((1:6)+6*(nodenum1-1)) ((1:6)+6*(nodenum2-1))],sk+(1:6))=conmat';
	  K(sk+(1:6),[((1:6)+6*(nodenum1-1)) ((1:6)+6*(nodenum2-1))])=conmat;
	  %We keep the first node in the model
	  %slavedofs=[slavedofs;6*(j-1)+(1:6)'];
	  slavedofs=[slavedofs;slave(slavedofs,[nodenum2 nodenum1],(1:6)')];
	  lines=[lines;nodenum1 nodenum2]; 
	end
	%slavedofs,pause
    else
      if length(findstr(b,'clamp'))==1;% Applying clamp conditions
	for j=1:numccs
	  sk=size(K,1);
	  nodenum1=data(j,1);
	  nodenum2=data(j,2);
	  conmat=[eye(6) -eye(6)];
	  K(sk+(1:6),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
	  K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:6))=conmat';
	  slavedofs=[slavedofs;((1:6)'+(nodenum2-1)*6)];%slavedofs,
							%we keep the
							%first coordinate
	end
      elseif length(findstr(b,'pin'))==1
	for j=1:numccs
	  sk=size(K,1);
	  nodenum1=data(j,1);
	  nodenum2=data(j,2);
	  udv=data(3:5);
	  udv=udv/norm(udv); % Normalize vector
	  [mv,mvl]=min(abs(udv)); % Find lowest component of vector
	  iv1=[0 0 0];iv1(mvl)=1; % Make this greatest in new
				  % independent vector
				  pervec1=cross(udv,iv1); % Cross product give perp to original vector
				  pervec2=cross(pervec1,udv); % And another perp to original vector
				  conmat=[eye(3) zeros(3)];
				  conmat=[conmat; zeros(2,3) [pervec1;pervec2]];
				  conmat=[conmat -conmat];ccnum
				  K(sk+(1:5),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
				  K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:5))=conmat';
				  slavedofs=[slavedofs;(1:3)'+(nodenum2-1)*6];
				  slavedofs=[slavedofs;(nodenum2-1)*6+3+dropindeces(pervec1,pervec2)];
	end
      elseif length(findstr(b,'roller'))==1
	for j=1:numccs
	  sk=size(K,1);
	  nodenum1=data(j,1);
	  nodenum2=data(j,2);
	  rudv=data(3:5); %Rotation Unit direction vector
	  tudv=data(6:8);
	  rudv=rudv/norm(rudv); % Normalize vector      
	  tudv=tudv/norm(tudv); % Normalize vector
	  if abs(dot(rudv,tudv))>2*eps
	    warndlg('Roller vectors not perpendicular. Fixing translation.','Bad Constraint')
	    perp=cross(rudv,tudv);
	    tudv=cross(perp,rudv);
	  end
	  perp=cross(rudv,tudv);
	  conmat=[[rudv;perp] zeros(2,3);...
		  zeros(2,3) [tudv perp]];
	  conmat=[conmat -conmat];
	  K(sk+(1:4),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
	  K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:4))=conmat';
	  slavedofs=[slavedofs;...
		     (nodenum2-1)*6+dropindices(perp,rudv);...
		     (nodenum2-1)*6+3+dropindices(perp,tudv)];
	end
      elseif length(findstr(b,'surfaceball'))==1
	for j=1:numccs
	  sk=size(K,1);
	  nodenum1=data(j,1);
	  nodenum2=data(j,2);
	  tudv=data(3:5);
	  tudv=tudv/norm(tudv);
	  conmat=[tudv 0 0 0];
	  conmat=[conmat -conmat];
	  K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:1))=conmat';
	  K(sk+(1:1),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
	  [val,valindex]=max(abs(tudv));
	  slavedofs=[slavedofs;(nodenum2-1)*6+valindex];
	end
      elseif length(findstr(b,'ball'))==1
	conmat=[eye(3) zeros(3)];
	for j=1:numccs
	  sk=size(K,1);
	  nodenum1=data(j,1);
	  nodenum2=data(j,2);
	  conmat=[eye(3) zeros(3) -eye(3) zeros(3)];
	  K(sk+(1:3),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
	  K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:3))=conmat';
	  slavedofs=[slavedofs;(nodenum2-1)*6+[1;2;3]];
	end
      elseif length(findstr(b,'rbeam'))==1
	for j=1:numccs
	  sk=size(K,1);
	  nodenum1=data(j,1);
	  nodenum2=data(j,2);
	  tudv=nodes(data(j,2),:)-nodes(data(j,1),:);
	  % tangent unit direction vector is a bad variable name, but it's not being used elsewhere ...
	  rx=tudv(1);
	  ry=tudv(2);
	  rz=tudv(3);
	  skewsym=[0   rz  -ry;
		   -rz 0    rx;
		   ry  -rx   0];
	  conmat=[eye(6,6) [-eye(3,3) skewsym;zeros(3,3) ...
			    -eye(3,3)]];
	  K([((1:6)+6*(nodenum1-1)) ((1:6)+6*(nodenum2-1))],sk+(1:6))=conmat';
	  K(sk+(1:6),[((1:6)+6*(nodenum1-1)) ((1:6)+6*(nodenum2-1))])=conmat;
	  %slavedofs=[slavedofs;(nodenum2-1)*6+(1:6)'];
	  slavedofs=[slavedofs;slave(slavedofs,[nodenum1 nodenum2],(1:6)')];
	  lines=[lines;nodenum1 nodenum2]; 
	end
      elseif length(findstr(b,'surface'))==1|length(findstr(b,'rod'))== 1
	if length(findstr(b,'rod'))==1
	  data(j,3:5)=nodes(data(j,1))-nodes(data(j,2));
	end
	for j=1:numccs
	  sk=size(K,1);
	  nodenum1=data(j,1);
	  nodenum2=data(j,2);
	  tudv=data(3:5);
	  tudv=tudv/norm(tudv);
	  conmat=[tudv 0 0 0];
	  conmat=[conmat;zeros(3) eye(3)];
	  conmat=[conmat -conmat];
	  K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:4))=conmat';
	  K(sk+(1:4),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
	  [val,valindex]=max(abs(tudv));
	  slavedofs=[slavedofs;...
		     (nodenum2-1)*6+valindex;...
		     (nodenum2-1)*6+3+[1;2;3]];
	  lines=[lines;nodenum1 nodenum2]; 
	end
      end
    end
  end
    % Apply boundary conditions to stiffness matrix now done. Sizes
    % of other matrices must be modified accordingly
  sizk=size(K,1);
  if size(M,1)<sizk
    M(sizk,sizk)=0;
    Fepsn(sizk)=0;
    Fnln(sizk)=0;
    F(sizk)=0;
    Bthermal(sizk,1)=0;
  end
elseif strcmp('nastran',a)
  if exist('slavedofs')==0
    slavedofs=[];
  end
  % This is where we apply all of the boundary conditions for NASTRAN. 
  for i=1:length(ccs)
    ccnumb=i;% This is the index number of the
             % boundary condition under consideration
    b=ccs(i).type;
    data=ccs(i).ccdata;
    numccs=size(data,1); %Number of boundary conditions
    if length(findstr(b,'clamp'))==1;% Applying clamp conditions
      for j=1:numccs
        sk=size(K,1);
        nodenum1=data(j,1);
        nodenum2=data(j,2);
        %conmat=[eye(6) -eye(6)];
        %K(sk+(1:6),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
        %K([(1:6)+6*(nodenum1-1)
        %(1:6)+6*(nodenum2-1)],sk+(1:6))=conmat';
	
	eid=length(element)+1;
	element(eid).nodes=[nodenum1, nodenum2];
	element(eid).elementtype='rbar';
        fprintf(fido,'RBAR, %g, %g, %g, 123456,,,\n',[eid, nodenum1, nodenum2]);
	
	slavedofs=[slavedofs;((1:6)'+(nodenum2-1)*6)];%slavedofs,
                                                      %we keep the
                                                      %first coordinate
      end
    elseif length(findstr(b,'pin'))==1
      disp('Warning, pin-beams not supported')
      for j=1:numccs
        sk=size(K,1);
        nodenum1=data(j,1);
        nodenum2=data(j,2);
        udv=data(3:5);
        udv=udv/norm(udv); % Normalize vector
        [mv,mvl]=min(abs(udv)); % Find lowest component of vector
        iv1=[0 0 0];iv1(mvl)=1; % Make this greatest in new
                                % independent vector
        pervec1=cross(udv,iv1); % Cross product give perp to original vector
        pervec2=cross(pervec1,udv); % And another perp to original vector
        conmat=[eye(3) zeros(3)];
        conmat=[conmat; zeros(2,3) [pervec1;pervec2]];
        conmat=[conmat -conmat];
        K(sk+(1:5),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
        K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:5))=conmat';
        slavedofs=[slavedofs;(1:3)'+(nodenum2-1)*6];
        slavedofs=[slavedofs;(nodenum2-1)*6+3+dropindeces(pervec1,pervec2)];
      end
    elseif length(findstr(b,'roller'))==1
      disp('warning, roller constrain not supported')
      for j=1:numccs
        sk=size(K,1);
        nodenum1=data(j,1);
        nodenum2=data(j,2);
        rudv=data(3:5); %Rotation Unit direction vector
        tudv=data(6:8);
        rudv=rudv/norm(rudv); % Normalize vector      
        tudv=tudv/norm(tudv); % Normalize vector
        if abs(dot(rudv,tudv))>2*eps
          warndlg('Roller vectors not perpendicular. Fixing translation.','Bad Constraint')
          perp=cross(rudv,tudv);
          tudv=cross(perp,rudv);
        end
        perp=cross(rudv,tudv);
        conmat=[[rudv;perp] zeros(2,3);...
                zeros(2,3) [tudv perp]];
        conmat=[conmat -conmat];
        K(sk+(1:4),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
        K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:4))=conmat';
        slavedofs=[slavedofs;...
                   (nodenum2-1)*6+dropindices(perp,rudv);...
                   (nodenum2-1)*6+3+dropindices(perp,tudv)];
      end
    elseif length(findstr(b,'ball'))==1
      disp('warning, ball constraint not supported')
      conmat=[eye(3) zeros(3)];
      for j=1:numccs
        sk=size(K,1);
        nodenum1=data(j,1);
        nodenum2=data(j,2);
        conmat=[eye(3) zeros(3) -eye(3) zeros(3)];
        K(sk+(1:3),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
        K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:3))=conmat';
        slavedofs=[slavedofs;(nodenum2-1)*6+[1;2;3]];
      end
    elseif length(findstr(b,'surfaceball'))==1
      disp('warning, ball constraint not supported')
      for j=1:numccs
        sk=size(K,1);
        nodenum1=data(j,1);
        nodenum2=data(j,2);
        tudv=data(3:5);
        tudv=tudv/norm(tudv);
        conmat=[tudv 0 0 0];
        conmat=[conmat -conmat];
        K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:1))=conmat';
        K(sk+(1:1),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
        [val,valindex]=max(abs(tudv));
        slavedofs=[slavedofs;(nodenum2-1)*6+valindex];
      end
    elseif length(findstr(b,'rbeam'))==1
      for j=1:numccs
        sk=size(K,1);
        nodenum1=data(j,1);
        nodenum2=data(j,2);
	eid=length(element)+1;
	element(eid).nodes=[nodenum1, nodenum2];
	element(eid).elementtype='rbar';
        fprintf(fido,'RBAR, %g, %g, %g, 123456,,,\n',[eid, nodenum1, nodenum2]);
	
	slavedofs=[slavedofs;((1:6)'+(nodenum2-1)*6)];%slavedofs,
                                                      %we keep the
                                                      %first coordinate
      end
    elseif length(findstr(b,'surface'))==1|length(findstr(b,'rod'))== 1
      disp('warning, constraint rod not supported')
      if length(findstr(b,'rod'))==1
        data(j,3:5)=nodes(data(j,1))-nodes(data(j,2));
      end
      for j=1:numccs
        sk=size(K,1);
        nodenum1=data(j,1);
        nodenum2=data(j,2);
        tudv=data(3:5);
        tudv=tudv/norm(tudv);
        conmat=[tudv 0 0 0];
        conmat=[conmat;zeros(3) eye(3)];
        conmat=[conmat -conmat];
        K([(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)],sk+(1:4))=conmat';
        K(sk+(1:4),[(1:6)+6*(nodenum1-1) (1:6)+6*(nodenum2-1)])=conmat;
        [val,valindex]=max(abs(tudv));
        slavedofs=[slavedofs;...
                   (nodenum2-1)*6+valindex;...
                   (nodenum2-1)*6+3+[1;2;3]];
        lines=[lines;nodenum1 nodenum2]; 
      end
    end
  end
  % Apply boundary conditions to stiffness matrix now done. Sizes
  % of other matrices must be modified accordingly
%   size(K)
%   size(M)
%   pause
end



function dropindnums=dropindeces(vector1,vector2)
% Optimally picks between x,y, and z to drop the most predominant
% two axes vectors in the two vectors.

[val,valind]=max(abs(cross(vector1,vector2))); % This gives us
                                               % the axis most
                                               % parallel to the
                                               % pin. We'll keep
                                               % that rotation
                                               % and ditch the
                                               % others.

dropindexes=abs(([1;2;3]==ones(3,1)*valind)-1).*[1;2;3];
dropindnums=sort(dropindexes);
dropindnums=dropindnums(2:3);


function newslaves=slave(slavedofs,nodes,dofs)
% Given potential nodes to be slaved, slave one or the other
% based on whether one has already been slaved. A node can't be
% slaved twice, or it will cause the stiffness matrix to end up
% under-reduced. If both nodes have been zeroed, we have a
% potential problem of conflicting constraints, but will warn and
% just send back a null vector
    
newslaves=(nodes(1)-1)*6+dofs;
if sum(ismember(slavedofs,newslaves))>0
  newslaves=(nodes(2)-1)*6+dofs;
end
if sum(ismember(slavedofs,newslaves))>0
  newslaves=[];
  uiwait(warndlg(['Both nodes ' num2str(nodes(1)) ' and ' num2str(nodes(2)) ...
           ' have previously been reduced out of the model by a preveous ' ...
           'constraint. I expect failure of the analysis. Please do not ' ...
           'connect a large number of constraintes to one another.'] ...
          ,'Bad constraints.','modal'));
end
