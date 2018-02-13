function outputfinal=wfem2nastran(fname,settings,param)
% WFEM2NASTRAN(FILENAME,SETTINGS,PARAMETER) 
% SETTINGS(1)=1: reload from previous run (default 0)
% SETTINGS(2)=1: restart from previous run (default 0)
% SETTINGS(3)=1: reload and run with same actions (default 1)
% A reload skips the generation of the mass and stiffness matrices,
% instead loading them from the .mat file associated with filename
% as long as the .mat file is newer than the originating file. 
% A restart also loads results from a previous run. A restart is a
% reload plus more.
% PARAMETER is user defined values accessible from within the
% input file. The variable PARAMETER can be used instead of
% variable definitions in the input file. For instance, Young's
% modulus could be defined as PARAMETER(1), the shear modulus as
% PARAMETER(2), and so on. This allows values to be defined during the
% WFEM call instead of being explicitly defined in the input
% file. 
  
% The design philosophy is that the central code
% sends generic commands to element subroutines to 
% what the element routines do, not knowing anything 
% about the element. Element names are the same as 
% the m-file names.

% The element routines do all processing. They:
%     a) Make the D and B matrices
%     b) Generate the local element matrix
%     c) Convert the local element matrix 
%     d) Into global coordinates
%     e) Insert the global elemental matrices into
%             the global mass, stiffness, and other matrices (incl forces)
%     f) Put drawing information into the global
%             drawing matrices
%     g) Obtain stresses and strains for the chosen element
%     h) Put stress and strain drawing information into
%             global drawing matrices
%
% Them main code, WFEM, or it's subroutines, do everything else,
% including numerical integration, eigensolutions, application of
% boundary conditions, determination of reaction forces and
% drawing. Exactly how to do this is illustrated in this in beam3 and
% other element files. No knowledge of the types of elements is
% required by wfem. 
% 
% Elements must have the following modes. They do not all have to
% do anything. They just have to exist to avoid errors:
% make: make mass and stiffness matrices, assemble into global
% makes: make geometric stiffness matrices, assemble into globel
% str: obtain element stresses and strains from deflections
% 
%
% This philosophy is a little unusual, but it allows for the
% introduction of new elements without editing of any existing 
% WFEM code. This prevents introduction of bugs by accident, 
% and makes the code easier to use in a classroom environment
% where students need to learn FEM coding without getting 
% lost in the details of reading files, writing files, graphics...
% Of course, the full source code can always be viewed and 
% modified. If you have any improvements that you think I should 
% include in future release, please send them on.
% Copyright Joseph C. Slater, 7/26/2002.
% joseph.slater@wright.edu





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Global variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear global
global restart % Are restarts allowed (from pevious run... to save time)?
global reload
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
global bcs
global ccs
global fs
global fms
global filename
global actions
global ismatnewer
global slavedofs
global parameter
global fido
%Boundary Conditions
if exist('param')
  parameter=param;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Local variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output=[];
%if nargin==0
%  fname='example.txt';
%end
if nargin==2|nargin==3
  if length(settings)==3
    reload=settings(1);
    restart=settings(2);
    keepactions=settings(3);    
  elseif length(settings)==2
    reload=settings(1);
    restart=settings(2);
    keepactions=1;
  elseif length(settings)==1
    reload=settings;
    restart=0;
    keepactions=1;
  else
    reload=0;
    restart=0;
    keepactions=1;
    if length(settings)~=0
      uiwait(warndlg('Input setting vector too long.'),'modal')
    end
  end	
else
  reload=0;
  restart=0;
end
reload=0;
restart=0;
exist(fname)
exist(fname)==2


if exist(fname)~=2
end

%if 1==1%nargin>0

oldpath=pwd;
if nargin>0
  if exist(fname)~=2
    disp(['Filename ' fname ' doesn''t exist in ' pwd '.'])
    fname=input('Enter filename: ','s')
  end
end

if exist(fname)~=2
  if exist('uigetfile')==5
    [fname,pathname]=uigetfile('*.txt;*.msh','Select FEM file.');
  else
    while exist(fname)~=2
      disp(['Filename ' fname ' doesn''t exist in ' pwd '.'])
      fname=input('Enter filename: ','s')
    end
  end
  filename=fname;
  pathtofile=which(filename);
%  (pathtofile(1:findstr(pathtofile,fname)-1));
  cd(pathtofile(1:findstr(pathtofile,fname)-1));
  pwd;
%else
%  [filename,pathname]=uigetfile('*.txt;*.msh','Select FEM file.');
%  oldpath=pwd;
%  cd(pathname);
end
filename=fname;






filenameroot=filename(1:findstr(filename,'.')-1);
% Call the subroutine that parses and loads the finite element model
% into the data structure.
%exist([filenameroot '.mat'])

% $$$ if exist([filenameroot '.mat'])==2&(restart==1|reload==1)
% $$$ %  filename;
% $$$ %  filenameroot;
% $$$   sourcefileinfo=dir(filename);
% $$$ %  sourcefileinfo.date;
% $$$   savefileinfo=dir([filenameroot '.mat']);
% $$$ %  disp('Reloading from .mat file.')
% $$$ else
% $$$   savefileinfo=dir(filename);
% $$$   sourcefileinfo=savefileinfo;
% $$$ end

%savefileinfo;
%sourcefileinfo.date,pause
%datenum(savefileinfo.date)
%datenum(sourcefileinfo.date)
%savefileinfo
%datenum(savefileinfo.date)
%ismatnewer=(datenum(savefileinfo.date)>datenum(sourcefileinfo.date)& ...
%      exist([filenameroot '.mat'])==2);



%Probably only need loadfem here
% $$$ if (restart==1|reload==1)&ismatnewer
% $$$   eval(['load ' filenameroot ' element K Ks M nodes beamprops ' ...
% $$$         'matprops elprops lines surfs bcs ccs F actions' ...
% $$$         ' slavedofs']);
% $$$   disp('Reloading from .mat file')
% $$$   if keepactions==0
% $$$     clear actions
% $$$   end
% $$$ else
  loadfem(1)
  % Make elemental stiffness matrices and assemble
%  assemblemk
  % Apply the boundary conditions
%  applybcs('apply')
  % Apply constraint conditions
%  applyconstraint('apply')
  % In case we want to restart.
%  eval(['save ' filenameroot]);
%end

%istrainforces


doactions=1;
ii=0;
doneactions=0;
actions;


while doactions==1
  if exist('actions')&(length(actions)>0)&doneactions==0
    ii=ii+1;
    currentaction=actions{ii};
    noactionwarning=1;
    doneactions=1;% Only one action can be put into the NASTRAN deck.
% $$$     if ii==length(actions)
% $$$       doneactions=1;
% $$$     end
  else
    if exist('noactionwarning')~=1
      disp(sprintf(['\n\n\nI presume you want to do SOMETHING other than just ' ...
               'generate some matrices.\n']))
      disp(sprintf(['May I suggest one of the following:\n modalanalysis'...
               '\n staticanalysis\n plotundeformed\n' ...
               'plotinitialstrain\n']))
      disp(sprintf(['These commands can be put at the end of your input ' ...
               'file following an action command\n']))
      disp(sprintf(['SOL # left blank\n']))
      %fprintf(fido,'SOL  $insert solution type to left of dollar sign')
      noactionwarning=1;
    end
    currentaction='end'; % We can move on.
    %currentaction=input('What action would you like to take (''end'' to end)?\n','s');
  end
  
  if strcmp(currentaction,'findinitialstrain')
    %Not enabled yet.
    X=initialstrain;
    initialstrainfound=1;
  elseif strcmp(currentaction,'plotinitialstrain')
    %Not enabled yet
    if exist('initialstrainfound')==1
      if initialstrainfound~=1
        X=initialstrain;
      end
    else
      X=initialstrain;
    end
    plotgeom('deformed',nodes,lines,X)
  elseif strcmp(currentaction,'plotundeformed')
    %Not enabled yet
    if exist('surfs')==0
      if exist('lines')==0
        plotgeom('undeformed',nodes)
      else
        plotgeom('undeformed',nodes,lines)
      end
      else
        plotgeom('undeformed',nodes,lines,surfs)
    end
  elseif strcmp(currentaction,'plotdeformed')
    %Not enabled yet
    if ~exist('X')
      disp('You need some analysis first.')
    else
      if exist('surfs')==0
        if exist('lines')==0
          plotgeom('deformed',nodes,X)
        else
          plotgeom('deformed',nodes,lines,X)
        end
      else
        plotgeom('deformed',nodes,lines,surfs,X)
      end
    end
% $$$   elseif strcmp(currentaction,'planefit')
% $$$     if ~exist('X')
% $$$       disp('You need some analysis first.')
% $$$     else
% $$$       [rmsvalue,a,b,c]=planefit(X);
% $$$     end
  elseif strcmp(currentaction,'modalanalysis');
    %Must be enabled yet &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if isempty(fs)& restart==1
      eval(['load ' filenameroot ' fs fms']);
    end
    if isempty(fs)
      [fs,fms]=modalanalysis;
      eval(['save ' filenameroot]);
    end
    output.fs=fs;
  elseif strcmp(currentaction,'modalreview')
    %Not enabled yet
    if restart==1
      eval(['load ' filenameroot ' fs fms'])
    end
    if isempty(fs)% maybe it tried to load the results, but they
                  % weren't in the file! 
      disp('You need to perform the action modalanalysis first.')
      %[fs,fms]=modalanalysis;
      %eval(['save ' filenameroot]);
    else
      modalreview;
    end
  elseif strcmp(currentaction,'reducedofs')
    %Not enabled yet
    [Tr,Mr,Kr]=reducedofs;
  elseif strcmp(currentaction,'planefit')
    %Not enabled yet
    if exist('X')~=1;
      X=staticanalysis;
    end
    [rmsvalue,a,b,c]=planefit(X);
    output.rmsvalue=rmsvalue;
  elseif strcmp(currentaction,'staticanalysis')
    %Not enabled yet
    X=staticanalysis;
  elseif strcmp('',currentaction)||isempty(currentaction)||strcmp('%',currentaction(1))
  % do nothing
  elseif strcmp(currentaction,'end')|strcmp(currentaction, ...
                                            '''end''')
    doactions=0;
  else
    currentaction
    disp(sprintf(['I don''t understand what ' currentaction [' means, but ' ...
                        'I''ll try to do it anyway.']]))
    if strcmp(currentaction,'exit')
      disp(sprintf(['I don''t think you want to exit Matlab, so I''m ' ...
               'ignoring that.\n']))
    else
      try
        if strcmp(currentaction,'')
          currentaction='%';
        end
        eval(currentaction)
      catch
        lasterr
      end
    end
  end
end
if nargout==1
  outputfinal=output;
end




%Open nastran file for writing
%%%%
%%%% Need to add test for existence of file.
[fido,msgo]=fopen([filenameroot '.dat'],'w','native');
%Write out header for post processing
fprintf(fido,['ASSIGN OUTPUT2=''' filenameroot '.op2'' UNIT=12,' ...
		    ' UNFORMATTED\n']);
%name for reference
fprintf(fido,['ID ' filenameroot '\n']);
%
%%%%%%%%%%%%%%%%
%Type of analysis
analysistype=3;
disp('What number for analysis type?')
fprintf(fido,['SOL, ' num2str(analysistype) '\n']);
%
fprintf(fido,['TIME, 10\n']);
%end header part 1
fprintf(fido,'CEND\n');
%
%Title for output printings
fprintf(fido,['TITLE=' filenameroot '\n']);
%
%
fprintf(fido,'METHOD=10\n');
fprintf(fido,'DISP=ALL\n');
%fprintf(fido,'\n');

fprintf(fido,'BEGIN BULK\n');

fprintf(fido,'PARAM, POST, -2$needed for Ideas post-processing\n');
numberofmodes=num2str(floor(min([size(nodes,1)*3,30])));
fprintf(fido,['EIGR, 10, LAN,0.0,1000.0,'  numberofmodes ',' numberofmodes '$eigensolution method\n']);
disp('Maximum solvable natural frequency is 1000 Hz.')
% Dump nodes to NASTRAN GRIDS
for ii=1:size(nodes,1)
  fprintf(fido,'GRID,%g,,%8.4E,%8.4E,%8.4E,,, \n',[ii,nodes(ii,1),nodes(ii,2),nodes(ii,3)]);
%  fprintf(fido,['GRID,' num2str(ii) ', ' num2str(nodes(ii,1)) ',' num2str(nodes(ii,2)) ',' num2str(nodes(ii,3)) '\n']);
end

% Dump points to NASTRAN GRIDS

%for ii=1:size(points,1)
%  fprintf(fido,'GRID,%g,,%E,%E,%E,,, \n',[ii+size(nodes,1),points(ii,1),points(ii,2),points(ii,3)]);
%  fprintf(fido,['GRID,' num2str(ii+size(nodes,1)) ', ' num2str(points(ii,1)) ',' num2str(points(ii,2)) ',' num2str(points(ii,3)) '\n']);
%end



usedprops=[];
for ii=1:length(element)
  %need to now define, for each element, a call that does a nastran dump
%  ii 
%  usedprops
  propnum=element(ii).properties;
  if (sum(usedprops==propnum))==0
    %need to call to write out prop type and info, then append to usedprops
    fcall=[element(ii).elementtype '(''nastranpropertydump'',' num2str(ii) ',' ...
	 num2str(fido) ')']
    eval(fcall);
    usedprops=[usedprops propnum];
  end
  fcall=[element(ii).elementtype '(''nastrandump'',' num2str(ii) ',' ...
	 num2str(fido) ')'];
  eval(fcall);
end

applyconstraint('nastran') %This has a looping through all
                           %constraints already. It created a
                           %unique element ID by adding onto
                           %existing element numbers (creating a
                           %dummy element for each constraint

fprintf(fido,'ENDDATA\n');
fclose(fido)
%%%%%%%%%%%

cd(oldpath)
