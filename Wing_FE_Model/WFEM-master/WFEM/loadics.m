function loadics
% This function generates initial conditions for a time simulation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Global variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global actions
global restart
global restartsamemk
global K
global filename
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
global X0 %Initial Displacement Conditions
global V0 %Initial Velocity Conditions
global ccs %Constraint conditions
global parameter
global curlineno
global DoverL
global parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Local variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenameroot=filename(1:findstr(filename,'.')-1);%pause


units;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
[fid,msg]=fopen(filename,'r','native');
% Assume at least one line exists in the file
line='0';
ii=1;
elprops=[];
eltypenum=0;
eltype=[blanks(20)];
line=fgetl(fid);
curlineno=curlineno+1;
nodenum=1;
pointnum=1;
mode=0;
kill=[];
elpropnum=0;
elnum=0;
points=[0 0 0];
nodes=[];
% In case we have functional initial condition, we will need these
% values. 
x=nodes(:,1);y=nodes(:,2);z=nodes(:,3);

while isstr(line)
  line;
    ii=ii+1;
    if strcmp(deblank(line),'')
      line='%';
    end
    if strcmp(line(1),'%')
		       %disp('Comment line or blank line.')
    else
      % Parse the data file here for all the information in it.
      if length(findstr(lower(line),'element'))==1&length(findstr(lower(line),'propert'))==1
        %disp('Reading Element Properties.')
        mode=1;
      elseif length(findstr(lower(line),'points'))==1
%        disp('Reading point locations')
        mode=2;
      elseif length(findstr(lower(line),'nodes'))==1
%        disp('Reading node locations')
        mode=3;
      elseif isempty(findstr(lower(line),'element'))~=1
%        elementtype=deblank(line(1:findstr(line,'element')-1));
        %disp(['Reading ' elementtype ' elements.'])
%        if exist(elementtype)~=2
%          disp(['Element type ' elementtype [' doesn''t exist. Please' ...
%                              ' make sure '] elementtype ['.m is in your' ...
%                              ' path']])
%          kill=1;
%          break
%        end
        mode=4;
%        eltypenum=eltypenum+1;
%        eltype(eltypenum,1:20)=[elementtype blanks(20- ...
%                                                   length(elementtype))];
      elseif length(findstr(lower(line),'variable'))==1
        disp('Reading variables.')
        mode=6;
      %
      elseif length(findstr(lower(line),'initial'))==1&length(findstr(lower(line),'conditions'))==1
        disp('Reading initial conditions.')
        mode=7;
      else
        %      line
        %      numline=eval(['[' line ']']);
        if mode==1
          %elpropnum=elpropnum+1;
          %elprops(elpropnum).a=numline;
        elseif mode==2
% $$$ 	points(pointnum,:)=numline(2:4);
% $$$ 	pointnum=pointnum+1;
        elseif mode==3
% $$$         numline,nodes,pause
	elseif mode==4%
        elseif mode==5%
        elseif mode==6%Read variables
          if exist(line(1:findstr(line,'=')-1))~=0
            uiwait(warndlg([ line(1:findstr(line,'=')-1) ' cannot be ' ...
                            'used as a variable.'],...
                           ['Don''t use variable ' line(1:findstr(line,'=')-1) ...
                           ],'modal'))
          else
            eval([line ';']);
          end
	elseif mode==7% Apply initial conditions
	  %This is where the real work takes place. We already have
          %defined x, y, and z, for each of the nodes, allowing
          %them to be used in fevals.
	  %x0 (lowercase x) is a temporary storage, each od
          %numbered column corresponds to an initial displacement,
          %each even to its velocity. Order: x y z tx ty tz
	  x0(:,[1 2])=eval(['[' line ']']);
	  dofnum=2;
	  while dofnum<7
	    line=fgetl(fid);
	    if strcmp(line(1),'%')
	    else
	      x0(:,[2*(dofnum-1)+[1 2]])=eval(['[' line ']']);
	      dofnum=dofnum+1;
	    end
	  end
	  %Assigning into single vector(s)
	  sm=size(M,1);
	  X0=zeros(sm,1);
	  V0=zeros(sm,1);
	  X0(1:6:(sm-5),1)=x0(:,1);
	  X0(2:6:(sm-4),1)=x0(:,3);
	  X0(3:6:(sm-3),1)=x0(:,5);
	  X0(4:6:(sm-2),1)=x0(:,7);
	  X0(5:6:(sm-1),1)=x0(:,9);
	  X0(6:6:(sm-0),1)=x0(:,11);
	  V0(1:6:(sm-5),1)=x0(:,2);
	  V0(2:6:(sm-4),1)=x0(:,4);
	  V0(3:6:(sm-3),1)=x0(:,6);
	  V0(4:6:(sm-2),1)=x0(:,8);
	  V0(5:6:(sm-1),1)=x0(:,10);
	  V0(6:6:(sm-0),1)=x0(:,12);
	  line=fgetl(fid);
	  while ~strcmp(deblank(line),'')
	    %We will read and interpret individual lines here
	    icdef=eval(['[' line ']']);
	    if length(ic==4)
	      X0((icdef(1)-1)*6+icdef(2),1)=icdef(3)+X0((icdef(1)-1)*6+icdef(2),1);
	      V0((icdef(1)-1)*6+icdef(2),1)=icdef(4)+V0((icdef(1)-1)*6+icdef(2),1);
	    else
	      if icdef(2)<7
		%Must be a displacement IC
		X0((icdef(1)-1)*6+icdef(2),1)=icdef(3)+X0((icdef(1)-1)*6+icdef(2),1);
	      else
		%Must be a displacement IC
		V0((icdef(1)-1)*6+icdef(2)-6,1)=icdef(3)+V0((icdef(1)-1)*6+icdef(2)-6,1);
	      end
	    end
	    line=fgetl(fid);
	  end
	end
      end
    end
    line=fgetl(fid);
    curlineno=curlineno+1;
  end
  fclose(fid);
  
  % Create zero force vectors (we want to make sure that they have
  % the right dimension in the rare case that we don't have any loads
  % of one type, or in the case of no boundary conditions (which
  % fixes the size).
  
 
  
