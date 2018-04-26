function loadfem
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
global ccs %Constraint conditions
global parameter
global curlineno
global DoverL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Local variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if nargin==0
%  filename='example.txt';
%end
%if exist(filename)==0
%  disp('Running default (example.txt).')
%  filename='example.txt';
  %disp(['Filename ' filename ' doesn''t exist in ' pwd '.'])
  %filename=input('Enter filename: ','s')
%end
%filename
filenameroot=filename(1:findstr(filename,'.')-1);%pause

% $$$ if exist([filenameroot '.mat'])==2
% $$$   sourcefileinfo=dir(filename);
% $$$   savefileinfo=dir([filenameroot '.mat']);
% $$$ else
% $$$   savefileinfo=dir(filename);
% $$$   sourcefileinfo=savefileinfo;
% $$$ end

%datenum(savefileinfo.date)>datenum(sourcefileinfo.date)
%exist([filenameroot '.mat'])==2
%restart==1


% $$$ if datenum(savefileinfo.date)>datenum(sourcefileinfo.date)& ...
% $$$       exist([filenameroot '.mat'])==2&restart==1
% $$$   eval(['load ' filenameroot '.mat'])
% $$$ else
units;
curlineno=0;
[fid,msg]=fopen(filename,'r','native');
line=fgetl(fid);
curlineno=curlineno+1;
ismeshed=1;
while isstr(line)&ismeshed==1
  if strcmp(line,'')
    line='%';
  end
  if ~strcmp(line(1),'%')
    if isempty(findstr(line,'element'))~=1 & isempty(findstr(line, ...
						  'bay'))~=1
      ismeshed=0;
    end
  end
  line=fgetl(fid);
  curlineno=curlineno+1;
end
fclose(fid);

if ismeshed==0
  disp('Meshing finite element model')
  filename=meshfem(filename);
  disp(['Meshed FEM is in file named' filename])
end
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
% Here we will read the file for the first time in order to know what all of the
% nodes are that exist. This is important because some elements (e.g.
% beam3) insert internal nodes if they are not defined explicitely. In
% picking the node number for the internal node, an error would occur if
% the automatically generated node number matched one that is defined
% later. The next two while loops are identical except that only node-read
% actions are taken in the first loop, and only other actions are taken
% in the later loop.


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
        disp('Reading point locations')
        mode=2;
      elseif length(findstr(lower(line),'nodes'))==1
        disp('Reading node locations')
        mode=3;
      elseif isempty(findstr(lower(line),'element'))~=1
        elementtype=deblank(line(1:findstr(line,'element')-1));
        %disp(['Reading ' elementtype ' elements.'])
        if exist(elementtype)~=2
          disp(['Element type ' elementtype [' doesn''t exist. Please' ...
                              ' make sure '] elementtype ['.m is in your' ...
                              ' path']])
          kill=1;
          break
        end
        mode=4;
        eltypenum=eltypenum+1;
        eltype(eltypenum,1:20)=[elementtype blanks(20- ...
                                                   length(elementtype))];
      elseif length(findstr(lower(line),'variable'))==1
        disp('Reading variables.')
        mode=6;
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
          try
	    numline=eval(['[' line ']']);  
	  catch
	    lasterr
	    disp(['Error on line ' num2str(curlineno) ' of input file:'])
	    line
	    disp(['My apologies, but I''m sure the line number given' ...
		  ' is wrong'])
	    pause
	  end
	  
	  if size(nodes,1)==numline(1)
            warndlg(['A previous element has already created a node ' ...
                     'number ' num2str(numline(1)) '. Please define ' ...
                     'this node earlier on the input ' ...
                     'file.'],'Node already exists.')
            nodes
            pause
            return
          else
            nodes(nodenum,:)=numline(2:4);
            nodenum=nodenum+1;
          end
        elseif mode==4
% $$$         elnum=elnum+1;
% $$$         element(elnum).elementtype=elementtype;
% $$$         fcall=[element(elnum).elementtype '(''generate'',[' line ...
% $$$                '],' num2str(elnum) ');'];
% $$$         eval(fcall);
% $$$         %eval(['eldata.' elementtype '=[' line '];']);
% $$$         %element(elnum).elementtype=elementtype;
% $$$         %element(elnum).nodes=numline;
        elseif mode==5
        elseif mode==6
          if exist(line(1:findstr(line,'=')-1))~=0
            uiwait(warndlg([ line(1:findstr(line,'=')-1) ' cannot be ' ...
                            'used as a variable.'],...
                           ['Don''t use variable ' line(1:findstr(line,'=')-1) ...
                           ],'modal'))
          else
            eval([line ';']);
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
  Fepsn=sparse(zeros(size(nodes,1)*6,1));
  F=Fepsn;
  Fnln=F;
 
  
  % Opening all over again so that we can read everything ``for the first time''.
  [fid,msg]=fopen(filename,'r','native');
  ii=0;
  line=fgetl(fid);
  curlineno=curlineno+1;
  % Here we read everything else.
  while isstr(line)
    ii=ii+1;
    if strcmp(deblank(line),'')
      line='%';
    end
    if strcmp(line(1),'%')
      %disp('Comment line or blank line.')
    else
      % Parse the data file here for all the information in it.
      if length(findstr(lower(line),'element'))==1&length(findstr(lower(line),'propert'))==1
        disp('Reading Element Properties.')
        mode=1;
      elseif strcmp(lower(line),'points')
        disp('Reading point locations')
        mode=2;
      elseif strcmp(lower(line),'nodes')
        %disp('Reading node locations')
        mode=3;
      elseif isempty(findstr(lower(line),'element'))~=1
        elementtype=deblank(line(1:findstr(line,'element')-1));
        disp(['Reading ' elementtype ' elements.'])
        %disp(['Reading elements for type ' elementtype '.'])
        if exist(elementtype)~=2
          disp(['Element type ' elementtype [' doesn''t exist. Please' ...
                              ' make sure '] elementtype ['.m is in your' ...
                              ' path']])
          kill=1;
          break
        end
        mode=4;
        eltypenum=eltypenum+1;
        eltype(eltypenum,1:20)=[elementtype blanks(20- ...
                                                   length(elementtype))];
      elseif length(findstr(lower(line),'fix'))==1
        %  Some form of fixed boundary conditions follow
        applybcs('read',line,fid)
      elseif length(findstr(lower(line),'constrain'))==1|length(findstr(lower(line),'constrian'))==1
        %  Two DOFsnodes are linked in some set or another of DOFs
        applyconstraint('read',line,fid)
      elseif length(findstr(lower(line),'load'))==1
        mode=5;
        disp('Reading applied loads.')
      elseif length(findstr(lower(line),'variable'))==1
        mode=6;
      elseif length(findstr(lower(line),'action'))==1
        mode=7; % We need to keep this in until we get rid of the
                % whole stinking mode idea.
        actionnumber=0;
        while isstr(line)&strcmp(deblank(line),'')~=1
          actionnumber=actionnumber+1;
          line=fgetl(fid);
	  curlineno=curlineno+1;	    
          if strcmp(line(1),'%')~=1&(line~=-1)
            actions{actionnumber}=line;
          end
        end
        
      else
        numline=str2num(line);
        if mode==1
          %elprops(size(elprops,1)+1,1:length(numline))=numline;
          elpropnum=elpropnum+1;
          line;
          while length(findstr(line,'...'))~=0
            dddloc=findstr(line,'...')-1;
            line=[line(1:dddloc(1)) ' ' fgetl(fid)];
	    curlineno=curlineno+1;	      
          end
          elementproperties=eval(['[' line '];']);
          elprops(elpropnum).a=elementproperties;
          
          %elprops(elpropnum).a=numline;
        elseif mode==2
          line;
          numline;
          points(pointnum,:)=numline(2:4);
          pointnum=pointnum+1;
        elseif mode==3
% $$$         numline,nodes,pause
% $$$ 	if size(nodes,1)==numline(1)
% $$$           warndlg(['A previous element has already created a node ' ...
% $$$                    'number ' num2str(numline(1)) '. Please define ' ...
% $$$                               'this node earlier on the input ' ...
% $$$                               'file.'],'Node already exists.')
% $$$           nodes
% $$$           pause
% $$$           return
% $$$         else
% $$$           nodes(nodenum,:)=numline(2:4);
% $$$           nodenum=nodenum+1;
% $$$         end
        elseif mode==4
          elnum=elnum+1;
          element(elnum).elementtype=elementtype;
          fcall=[element(elnum).elementtype '(''generate'',[' line ...
                 '],' num2str(elnum) ');'];
          eval(fcall);
          %eval(['eldata.' elementtype '=[' line '];']);
          %element(elnum).elementtype=elementtype;
          %element(elnum).nodes=numline;
        elseif mode==5
          nodenum=numline(1);
          dofnum=numline(2);
          localload=numline(3);
          F(dofnum+6*(nodenum-1))=localload;
        elseif mode==6
        end
      end
    end
    line=fgetl(fid);
    curlineno=curlineno+1;
    %,pause
  end
  fclose(fid);
%  eval(['save ' filenameroot '.mat']);
%end
%size(Fepsn)
%if exist('restart')==0
%  restart=0;
%end