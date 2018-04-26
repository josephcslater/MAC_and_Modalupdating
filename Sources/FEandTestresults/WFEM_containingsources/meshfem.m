function filenameout=meshfem(filenamein)

% Copyright Joseph C. Slater, 7/26/2002.
% joseph.slater@wright.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Global variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global K
global Ks
global M
%global bay
%global nodes
global beamprops
global matprops
global points
global Fepsn %Forces representing initial strain
global parameter
nodes=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Local variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
units
%if nargin==0
%  filenamein='compactexample.txt';
%end
filenameout=[filenamein(1:(findstr(filenamein,'.')-1)) '.msh'];
if strcmp(filenamein,filenameout)==1
  filenameout=[filenamein(1:(findstr(filenamein,'.')-1)) '.mes'];
  warndlg('Bad extension.',['Input file to be meshed should not have' ...
		    ' extension ''.msh''. ''.mes'' being used for' ...
		    ' meshed file. Please change extensions after meshing.'])
end
if exist(filenamein)==0
  disp(['Filename ' filenamein ' doesn''t exist in ' pwd '.'])
  filenamein=input('Enter filename: ','s')
end
%filenameout='out.txt'
[fidi,msg]=fopen(filenamein,'r','native');
[fido,msgo]=fopen(filenameout,'w','native');
% Assume at least one line exists in the file
line='0';
ii=1;
elprops=[];
eltypenum=0;
eltype=[blanks(20)];
line=fgetl(fidi)
nodenum=1;
pointnum=1;
mode=0;
kill=[];
while isstr(line)==1

  ii=ii+1;
  if strcmp(deblank(line),'')==1
    % Blank lines, just rewrite them out.
    line='\n';
    fprintf(fido,[line '\n']);
  elseif strcmp(lower(line),'nodes')==1
    % Nodes. We need to load these into memory in case we need
    % them, and we also need to write them back out. They must be
    % evaluated as well (Matlab isn't going to do symbolic math in
    % this code!)
    disp('Reading Nodes')
    isnode=1;
    fprintf(fido,'nodes\n');
    line=fgetl(fidi);
    while isnode==1
      % Read all of the nodes in the clump of nodes.
      ii=ii+1;%if length(line)==0, line='%';end
      %disp('node')
      %line
      %str2numline=str2num(line),
      %pause
      line;
      if length(deblank(line))>0&~strcmp(line(1),'%')
        line;
        isstr(line);
	try
	  str2numline=eval(['[' line '];']);
	catch
	  lasterr
	  sprintf(['This variable needs to be defined in your input '...
		 'file.\n You might as well Control-c now.'])
	  pause
	end
	
        %pause
	nodes=[nodes;str2numline(1,2:4)];
	fprintf(fido,[line '\n']);
      elseif strcmp(line(1),'%')
	fprintf(fido,['%' line '\n']);
      else
	isnode=0;
      end
      line=fgetl(fidi);
      if strcmp(line,'')
        fprintf(fido,'\n');
        isnode=0;
      end
    end

  elseif length(findstr(lower(line),'variable'))==1
    % We need to read in and evaluate parameters, in case they are
    % used. We write these lines out for reference only.
    fprintf(fido,['variables\n']);
    line=fgetl(fidi);
    fprintf(fido,[line '\n']);
    disp('Reading parameters')
    while length(deblank(line))>0
      %line,pause
      if exist(line(1:findstr(line,'=')-1))~=0
	%alpha
        uiwait(warndlg([ line(1:findstr(line,'=')-1) ' cannot be ' ...
                        'used as a variable.'],...
                       ['Don''t use variable ' line(1:findstr(line,'=')-1) ...
                       ],'modal'))
      else
        eval([line ';']);
      end
      line=fgetl(fidi);
      fprintf(fido,[line '\n']);
    end
  elseif length(findstr(lower(line),'action'))==1
    fprintf(fido,[line '\n']);
    line=fgetl(fidi);
    while isstr(line)
      fprintf(fido,[line '\n']);
      line=fgetl(fidi);
    end
  elseif length(findstr(lower(line),'global'))==1
%     % We need to read in and make globals available in case they
%     % are also used (like we did with variables). They are also
%     % written out only for reference. 
%     fprintf(fido,['globals\n']);
%     line=fgetl(fidi);
%     fprintf(fido,[line '\n']);
%     disp('Making globals')
%     while length(deblank(line))>0
%       %line,pause
%       if exist(deblank(line))
%         disp(['Variable ' line(1:findstr(line,'=')-1) ' cannot be ' ...
%               'used as a global variable.'])
%         beep,beep,beep
%         disp(['Please cancel this run with Control-C and fix the ' ...
%               'input file to use another global variable name.'])
%         pause
%       else
%         eval(['global' line ';']);
%       end
%       line=fgetl(fidi);
%       fprintf(fido,[line '\n']);
%     end
  elseif isempty(findstr(line,'element'))~=1 & isempty(findstr(line, ...
                                                      'bay'))~=1
    % This is the cool part. When we find a bay of elements we
    % process and repeat them as necessary to generate the meshed model.
    disp('Bay of elements found.')
    isbay=1;elementtype=0;eltypenum=0;
    clear bay
    while isbay==1
      % Get all of the bay information. Check that this makes
      % sense. 
      ii=ii+1;%pause
      line=fgetl(fidi);
      if strcmp(line,'')
        line='%';
      end
      if isempty(findstr(line,'element'))==1&elementtype==0
        if ~strcmp('%',line(1))
          warndlg(['Missing element definition. Line ' num2str(ii)],'Bad bay.')
        end 
      elseif findstr(line,'repeat')
        disp('Creating Bays.')
        isbay=0;
      elseif isempty(findstr(line,'element'))~=1
        eltypenum=eltypenum+1;elnum=0;
        elementtype=deblank(line(1:findstr(line,'element')-1));
        disp(['Reading, in block, elements of type ' elementtype '.'])
        if exist(elementtype)~=2
          warndlg(['Please make sure ' elementtype '.m is in your ' ...
                   'path'],['Element ' elementtype ' doesn''t exist.'])
          kill=1;
          break
        else
          bay(1).eltype(eltypenum).element=elementtype;%bay(1).eltype(eltypenum).element,pause
        end
      elseif strcmp(line(1),'%')
      else
        elnum=elnum+1;
        eldata= eval(['[' line ']']);%eval([line])  ;
        bay(1).eltype(eltypenum).properties(elnum,:)=eldata;
      end
    end
%    global bay
    %Parse the repeat line.
    line=lower(line);
    p1=findstr(line,'bay');p2=findstr(line,'times');
    p3=findstr(line,'ttach');p4=findstr(line,'to');
    p5=max([findstr(line,'.') ]);
%    line,p1,p2,p3,p4,p5
    %numbays=str2num(line((p1+3):(p2-1))); 
    %line((p1+3):(p2-1)),n,pause
    numbays=eval(['[' line((p1+3):(p2-1)) '];']);
    attachnodes1=str2num(line((p3+6):(p4-1)));% These are the start of each bay.
    attachnodes2=str2num(line((p4+2):(p5-1)));% These are free unless
                                    % connected to by the next bay.
    lanodes=length(attachnodes1);
    % Here we need to figure out, with the help of the element
    % subroutines, what actual nodes exist in this block
    baynodes=[];
    %pause
    for elnum=1:length(bay(1).eltype)
      fcall=[bay(1).eltype(elnum).element '(''numofnodes'',[' ...
                       num2str(bay(1).eltype(elnum).properties(1,:)) '])'];
      numofnodes=eval(fcall);
      elnodes=bay(1).eltype(elnum).properties(:,1:numofnodes);
      %disp('start')
      %baynodes
      %unique(elnodes)
      baynodes=[baynodes;unique(elnodes)];  
    end
    baynodes=unique(baynodes)';lbaynodes=length(baynodes);
%    baynodes,attachnodes1, attachnodes2
    internalnodes=setdiff(baynodes,[attachnodes1,attachnodes2]);
    %Need to create a new group of nodes, and set up the mapping.
    %attachnodes2,attachnodes1,nodes
    %pause

    % Here we check that all attachment nodes are displaced the
    % same amount. If not, this is an indicator that the attachment
    % wasn't defined properly.
    posdiffsvd=sort(-svd(nodes(attachnodes2,:)-nodes(attachnodes1,:)));
    
    if abs(posdiffsvd(2))>(1e-14)*max(max(nodes))
      warndlg(['Bad connectors in bay 2. Other bays not even checked.'],'Warning')
    end
    posdiff=nodes(attachnodes2(1),:)-nodes(attachnodes1(1),:);%pause% 

    if length(attachnodes1)~=length(attachnodes2)
      warndlg(['Attachments must have equal number of points. Line ' num2str(num2str(ii))],['Bad ' ...
                          'bay.'])
    else
      curmaxnodenum=size(nodes,1);
      locate(attachnodes1)=1:lanodes;
      locate(attachnodes2)=(1:lanodes)+lanodes;
      locate(internalnodes)=(lanodes*2+1):length(baynodes);
      baynodes=[attachnodes1 attachnodes2 internalnodes];
      for baynum=2:numbays+1%numbays is the number of repeat bays.
        disp(['Creating bay ' num2str(baynum) '.'])
        curmaxnodenum=size(nodes,1);
        baynodes(baynum,:)=[baynodes(baynum-1,(lanodes+1):lanodes*2) curmaxnodenum+(1:(lbaynodes-lanodes))];
        baynodes;%pause
        %baynodes
	[posdiffmesh,newnodes]=meshgrid(posdiff, curmaxnodenum+(1: ...
                                                          (lbaynodes-lanodes)));
        %nodes
        nodes( curmaxnodenum+(1:lbaynodes-lanodes),:)=posdiffmesh+...
            nodes( baynodes(baynum-1,lanodes+1:lbaynodes),:);% Create new nodes at correct locations.
% $$$         new numberings are already created here. Need to put them ...
% $$$             in the structured array positions.
        %nodes,pause
	fprintf(fido,['\n%%Bay ' num2str(baynum) ' nodes.\n']);
	fprintf(fido,['nodes\n']);
	for baynodenum=curmaxnodenum+(1:lbaynodes-lanodes)
	  fprintf(fido,[num2str(baynodenum) ' ' num2str(nodes( baynodenum,:)) '\n']);
	end
        %nodes and baynodes are fine 11/26/02 JCS
        bay(baynum)=bay(1);
	for eltypenum=1:length(bay(1).eltype)
	  %          bay(baynum)=bay(1);
	  
	  
% $$$ 	  [bay(1).eltype(eltypenum).element '(''numofnodes'',[' ...
% $$$ 		 num2str(bay(1).eltype(eltypenum).properties(1,:)) '])']
	  feval=[bay(1).eltype(eltypenum).element '(''numofnodes'',[' ...
		 num2str(bay(1).eltype(eltypenum).properties(1,:)) '])'];
          
          numofnodes=eval(feval);% we need to know how many,
                                 % starting from the front, of the
                                 % properties values are actually
                                 % node numbers.
          bn=baynodes(baynum,:);
          %eltypenum,pause
          %	  locate
%	  bay(baynum).eltype(elnum).properties;%(:,1:numofnodes)
%	  bn(locate(bay(baynum).eltype(elnum).properties(:,1:numofnodes)));%pause

bay(baynum).eltype(eltypenum).properties(:,1:numofnodes)=bn(locate(bay(baynum).eltype(eltypenum).properties(:,1:numofnodes)));
%	  bay(baynum).eltype(1).properties;%pause
%          bay(baynum).eltype(1).element
%          bay(baynum).eltype(1).properties
 %         bay(baynum).eltype(2).element
 %         bay(baynum).eltype(2).properties
        end
      end
%bay.eltype,pause
      % Here's where we start to write out the bay repeats, over,
      % an over...
      if length(bay)==1
        disp(['Livin'' on the edge with zero repeats, are we? No ' ...
              'problem. I work for nothing anyway.']) 
      end
%bay(1).eltype(1).element
%pause
      for baynum=1:length(bay)
	fprintf(fido,['\n%%Begin bay number ' num2str(baynum) '.\n']);

	for eltype=1:length(bay(baynum).eltype)
	  fprintf(fido,[bay(baynum).eltype(eltype).element ' elements\n']);
	  for elnum=1:size(bay(baynum).eltype(eltype).properties,1)
	    fprintf(fido, ...
		    [num2str(bay(baynum).eltype(eltype).properties(elnum,:)) '\n']);
	  end
          fprintf(fido,'\n');
	end
	fprintf(fido,['%%End bay number ' num2str(baynum) '.\n\n']);
      end
  
    end
  else
    if strcmp(line(1),'%')
      line=['%' line]; % I have no good reason for adding the extra
             % '%', but I do anyway.
    end
  % Write out the current line to the output file because it is
  % boring. This is what happens to anything else. Anything that
  % isn't picked up above. 
    fprintf(fido,[line '\n']);
  end      
  % Start all over with a new line, go back to the top, and try to
  % figure out what ramifications the line has.
  line=fgetl(fidi);
  
end
% Can't leave files open. It's a bad thing. 
fclose(fidi);
fclose(fido);
