function applybcs(a,b,c)
% APPLYBCS applies boundary conditions using the Lagrange
% Multiplier technique, see Cook, Malkus and Plesha. See wfem.pdf
% for details. 

%Copyright Joseph C. Slater, 2002
  
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
global bcs   % Boundary Conditions
global slavedofs
global Bthermal

%disp('hello')
if strcmp('read',a)
  % This is going to read in all of the boundary condition data.
  command=b;
  fid=c;
  line=fgetl(fid);
  while strcmp('%',line(1))
      line=fgetl(fid);
  end
  numberofnumbers=length(str2num(line));
  i=0;
  % We're loading all of the data in at once. We'll end up with
  % chunks of like boundary conditions, separated only if the user
  % does so. I'd like to separate all of these later if
  % possible. This will allow spiffier plotting.
  while length(str2num(line))==numberofnumbers
    i=i+1;line;
    bcdata(i,:)=str2num(line);
    line=fgetl(fid);line=[line ' '];
    while strcmp('%',line(1))
        line=fgetl(fid);line=[line ' '];
    end
  end
  
  if exist('bcs')
    bcnum=length(bcs)+1;
  else
    bcnum=1;
  end
  
  if length(findstr(b,'clamp'))==1
    bcs(bcnum).type='clamp';
    bcs(bcnum).bcdata=bcdata;
  elseif length(findstr(b,'pin'))==1
    bcs(bcnum).type='pin';
    bcs(bcnum).bcdata=bcdata;
  elseif length(findstr(b,'roller'))==1
    bcs(bcnum).type='roller';
    bcs(bcnum).bcdata=bcdata;
  elseif length(strfind(b,'surfaceball'))==1
    bcs(bcnum).type='surfaceball';
    bcs(bcnum).bcdata=bcdata;
  elseif length(findstr(b,'ball'))==1
    bcs(bcnum).type='ball';
    bcs(bcnum).bcdata=bcdata;
  elseif length(findstr(b,'surface'))==1
    bcs(bcnum).type='surface';
    bcs(bcnum).bcdata=bcdata;
  end

  
elseif strcmp('apply',a)
  if exist('slavedofs')==0
    slavedofs=[];
  end
  % This is where we apply all of the boundary conditions. 
  for i=1:length(bcs)
    bcnum=i;% This is the index number of the
             % boundary condition under consideration
    b=bcs(i).type;
	%b
    data=bcs(i).bcdata;
    numbcs=size(data,1); %Number of boundary conditions
    if length(findstr(b,'clamp'))==1;% Applying clamped boundary conditions
      for j=1:numbcs
        sk=size(K,1);
        nodenum=data(j,1);
        K((1:6)+6*(nodenum-1),sk+(1:6))=eye(6);
        K(sk+(1:6),(1:6)+6*(nodenum-1))=eye(6);
        slavedofs=[slavedofs;((1:6)'+(nodenum-1)*6)];%slavedofs,pause
      end
    elseif length(findstr(b,'pin'))==1
      for j=1:numbcs
        sk=size(K,1);
        nodenum=data(j,1);
	
        udv=data(j,2:4);
        udv=udv/norm(udv); % Normalize vector
        [mv,mvl]=min(abs(udv)); % Find lowest component of vector
        iv1=[0 0 0];iv1(mvl)=1; % Make this greatest in new
                                % independent vector
        pervec1=cross(udv,iv1); % Cross product give perp to original vector
        pervec2=cross(pervec1,udv); % And another perp to original vector
        conmat=[eye(3) zeros(3)];
        conmat=[conmat; zeros(2,3) [pervec1;pervec2]];
        K((1:6)+6*(nodenum-1),sk+(1:5))=conmat';
        K(sk+(1:5),(1:6)+6*(nodenum-1))=conmat;
        slavedofs=[slavedofs;(1:3)'+(nodenum-1)*6];
        slavedofs=[slavedofs;(nodenum-1)*6+3+dropindeces(pervec1,pervec2)];
      end
    elseif length(findstr(b,'roller'))==1
      for j=1:numbcs
        sk=size(K,1);
        nodenum=data(j,1);
        rudv=data(j,2:4); %Rotation Unit direction vector
        tudv=data(j,5:7);
        rudv=rudv/norm(rudv); % Normalize vector      
        tudv=tudv/norm(tudv); % Normalize vector
        if abs(dot(rudv,tudv))>2*eps
          warndlg('Roller vectors not perpendicular. Fixing translation.','Bad Constraint')
          perp=cross(rudv,tudv);
          tudv=cross(perp,rudv);
        end
        perp=cross(rudv,tudv);
        conmat=[[rudv;perp] zeros(2,3);...
                zeros(2,3) [tudv ;perp]];
        K((1:6)+6*(nodenum-1),sk+(1:4))=conmat';
        K(sk+(1:4),(1:6)+6*(nodenum-1))=conmat;
        slavedofs=[slavedofs;...
                   (nodenum-1)*6+dropindeces(perp,rudv);...
                   (nodenum-1)*6+3+dropindeces(perp,tudv)];
      end
    elseif length(strfind(b,'surfaceball'))==1
      for j=1:numbcs
        sk=size(K,1);
        nodenum=data(j,1);
		%disp('stuff now')
		%data
        tudv=data(j,2:4);
        tudv=tudv/norm(tudv);
        conmat=[tudv 0 0 0];
        K((1:6)+6*(nodenum-1),sk+(1:1))=conmat';
        K(sk+(1:1),(1:6)+6*(nodenum-1))=conmat;
        [val,valindex]=max(abs(tudv));
        slavedofs=[slavedofs;(nodenum-1)*6+valindex];
      end
    elseif length(findstr(b,'ball'))==1
      conmat=[eye(3) zeros(3)];
      for j=1:numbcs
        sk=size(K,1);
        nodenum=data(j,1);
        conmat=[eye(3) zeros(3)];
        K((1:6)+6*(nodenum-1),sk+(1:3))=conmat';
        K(sk+(1:3),(1:6)+6*(nodenum-1))=conmat;
        slavedofs=[slavedofs;(nodenum-1)*6+[1;2;3]];
      end
    elseif length(findstr(b,'surface'))==1
      for j=1:numbcs
        sk=size(K,1);
        nodenum=data(j,1);
        tudv=data(j,2:4);
        tudv=tudv/norm(tudv);
        conmat=[tudv 0 0 0];
        conmat=[conmat;zeros(3) eye(3)];
        K((1:6)+6*(nodenum-1),sk+(1:4))=conmat';
        K(sk+(1:4),(1:6)+6*(nodenum-1))=conmat;
        [val,valindex]=max(abs(tudv));
        slavedofs=[slavedofs;...
                   (nodenum-1)*6+valindex;...
                   (nodenum-1)*6+3+[1;2;3]];
       
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
                                                 % others
  
  dropindexes=abs(([1;2;3]==ones(3,1)*valind)-1).*[1;2;3];
  dropindnums=sort(dropindexes);
  dropindnums=dropindnums(2:3);
