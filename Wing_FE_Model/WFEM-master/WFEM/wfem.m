function outputfinal=wfem(fname,settings,param)
% WFEM(FILENAME,SETTINGS,PARAMETER) 
% FILENAME is a string containing the name of the input file.
% All other input variables are optional. 
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
global X0 V0
global filename
global actions
global ismatnewer
global slavedofs
global parameter
global Tr  Mr Kr
global Cd Cv Ca Bso Xd
global nin nout
global Bthermal

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
if nargin==2||nargin==3
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

if nargin>0
  if ~exist(fname)
    disp(['Filename ' fname ' doesn''t exist in ' pwd '.'])
    fname=input('Enter filename: ','s')
  end
  oldpath=pwd;
  if ~exist(fname)
    [fname,pathname]=uigetfile('*.txt;*.msh','Select FEM file.');
    cd(pathname)
  end
  filename=fname;
  pathtofile=which(filename);
else
  [filename,pathname]=uigetfile('*.txt;*.msh','Select FEM file.');
  oldpath=pwd;
  cd(pathname);
end

filenameroot=filename(1:findstr(filename,'.')-1);
% Call the subroutine that parses and loads the finite element model
% into the data structure.
%exist([filenameroot '.mat'])

if exist([filenameroot '.mat'])==2&(restart==1|reload==1)
%  filename;
%  filenameroot;
  sourcefileinfo=dir(filename);
%  sourcefileinfo.date;
  savefileinfo=dir([filenameroot '.mat']);
%  disp('Reloading from .mat file.')
ismatnewer=(datenum(savefileinfo.date)>datenum(sourcefileinfo.date)& ...
      exist([filenameroot '.mat'])==2);
else
  savefileinfo=dir(filename);
  sourcefileinfo=savefileinfo;
  ismatnewer=0;
end

%savefileinfo;
%sourcefileinfo.date,pause
%datenum(savefileinfo.date)
%datenum(sourcefileinfo.date)
%savefileinfo
%datenum(savefileinfo.date)
savefileinfo.date
%ismatnewer=(datenum(savefileinfo.date)>datenum(sourcefileinfo.date)& ...
%      exist([filenameroot '.mat'])==2);

if (restart==1|reload==1)&ismatnewer
  eval(['load ' filenameroot ' element K Ks M nodes beamprops ' ...
        'matprops elprops lines surfs bcs ccs F actions' ...
        ' slavedofs']);
  disp('Reloading from .mat file')
  if keepactions==0
    clear actions
  end
else
  loadfem
  %loadics
  disp('Load ICS disabled')
  % Make elemental stiffness matrices and assemble
  assemblemk
  % Apply the boundary conditions
  applybcs('apply')
  % Apply constraint conditions
  applyconstraint('apply')
  % In case we want to restart.
  eval(['save ' filenameroot]);
end

istrainforces;

%ndofs=size(nodes,1)*6;global ndofs;
%dofs=(1:ndofs)';global dofs;
%remove([1 2 3 4 5 7 8 9 10 11 13 15 16 17])
%eig(M)
%pause
%sqrt(eig(M\K))
%whos
%clf
%eig(K)
%modalanalysis
%staticanalysis
doactions=1;
ii=0;
doneactions=0;
actions;
while doactions==1
  if exist('actions')&(length(actions)>0)&doneactions==0
    ii=ii+1;
    currentaction=actions{ii};
    noactionwarning=1;
    if ii==length(actions)
      doneactions=1;
    end
  else
    if exist('noactionwarning')~=1
      disp(sprintf(['\n\n\nI presume you want to do SOMETHING other than just ' ...
               'generate some matrices.\n']))
      disp(sprintf(['May I suggest one of the following:\nmodalanalysis'...
               '\nstaticanalysis\nplotundeformed\' ...
               'nplotinitialstrain\n']))
      disp(sprintf(['These commands can be put at the end of your input ' ...
               'file following an action command\n']))
      noactionwarning=1;
    end
    currentaction=input('What action would you like to take (''end'' to end, see manual for others)?\nWFEM>>>  ','s');
  end
  
  if strcmp(currentaction,'findinitialstrain')
    X=initialstrain;
    initialstrainfound=1;
  elseif strcmp(currentaction,'plotinitialstrain')
    if exist('initialstrainfound')==1
      if initialstrainfound~=1
        X=initialstrain;
      end
    else
      X=initialstrain;
    end
    plotgeom('deformed',nodes,lines,surfs,X)
  elseif strcmp(currentaction,'plotundeformed')
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
    if isempty(fs)& restart==1
      eval(['load ' filenameroot ' fs fms']);
    end
    if isempty(fs)
      [fs,fms]=modalanalysis;
      eval(['save ' filenameroot]);
    end
    output.fs=fs;
  elseif strcmp(currentaction,'modalreview')
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
    [Tr,Mr,Kr]=reducedofs;
  elseif strcmp(currentaction,'planefit')
    if exist('X')~=1;
      X=staticanalysis;
    end
    [rmsvalue,a,b,c]=planefit(X);
    output.rmsvalue=rmsvalue;
  elseif strcmp(currentaction,'staticanalysis')
    X=staticanalysis;
  elseif strcmp('',currentaction)||isempty(currentaction)||strcmp('%',currentaction(1))
  % do nothing
  elseif strcmp(currentaction,'end')|strcmp(currentaction, ...
                                            '''end''')
    doactions=0;
  else
    currentaction;
    %disp(sprintf(['I don''t understand what ' currentaction [' means, but ' ...
    %                    '\n I''ll try to do it anyway.']]))
	%    fprintf('I don''t understand what {%s} means, but\nI''ll try to do it anyway.\n' ,currentaction  )
	if length(currentaction)>5&& strcmp(currentaction(1:4),'disp')
		else
			fprintf('   Trying: ''%s''.\n' ,currentaction  )
	end
	
			
	
	
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
cd(oldpath)
