function plotgeom(mode,nodes,lines,surfs,displacements,scale,T)
%Copyright Joseph C. Slater, 9/30/2002
% scale is a vector containing the scaling amplitude ([0,1]) and
% the shade value for plotting the undeformed wire mesh (default 0.75)
% These are configuration values to attempt to approximate
% reasonable maximum amplitudes of displacements for plotting.

  
maxrotplot=0.3;
maxdisplot=0.1;
maxsizex=max(nodes(:,1));
minsizex=min(nodes(:,1));
avex=(maxsizex+minsizex)/2;
maxsizey=max(nodes(:,2));
minsizey=min(nodes(:,2));
avey=(maxsizey+minsizey)/2;
maxsizez=max(nodes(:,3));
minsizez=min(nodes(:,3));
avez=(maxsizez+minsizez)/2;
maxnodespread=max([(maxsizex-minsizex),(maxsizey-minsizey), ...
                   (maxsizez-minsizez)]);

numnodes=size(nodes,1);
numlines=size(lines,1);
%displacements=displacements/max(abs(displacements));
if strcmp(mode,'undeformed')
  numnodes=size(nodes,1);
  if exist('lines')
    numlines=size(lines,1);
  else
    numlines=0;
  end
  
  if exist('surfs')
    numsurfs=size(surfs,1);
  else 
    numsurfs=0;
  end
  %for nodenum=1:numnodes
  h=plot3(nodes(:,1),nodes(:,2),nodes(:,3),'o');
  set(h,'markersize',3)
  cameratoolbar('show')
  hold on
  for nodenum=1:size(nodes,1)
    nodeh=text(nodes(nodenum,1),nodes(nodenum,2),nodes(nodenum,3), ...
               num2str(nodenum));
    set(nodeh,'color',[1 0 0],'fontsize',11,'margin',5)
%    set(nodeh,'buttondownfcn',['disp(get(' num2str(nodeh,16)
%    ',''string''))'])
    set(nodeh,'buttondownfcn',['set(' num2str(nodeh,16) ',''fontsize'',' ...
		    '(25-(get(' num2str(nodeh,16) ',''fontsize'' )-11) ))'])
    %nodeh
    %get(nodeh)
    %pause
  end
  
  for linenum=1:numlines
    curline=[nodes(lines(linenum,1),:);nodes(lines(linenum,2),:)];
    h=plot3(curline(:,1),curline(:,2),curline(:,3));
    %alpha(h,.7);
  end
  for surfnum=1:numsurfs
    coords=nodes(surfs(surfnum,1:4),:);
    patchcolor=surfs(surfnum,5:7);
    h=patch(coords(:,1),coords(:,2),coords(:,3),patchcolor);
    alpha(h,.7)
  end
  hold off
  grid on
  dnodes=nodes;
  
end
if exist('scale')==0
  scale=1;
end

if strcmp(mode,'deformed')
  %disp('plotdeformed')
  %nargin,pause
  if nargin==4
    displacements=surfs;
    clear surfs
    scale=1;
  elseif nargin==5
    if max(size(displacements))<(size(nodes,1)*6-1)
      scale=displacements;
      displacements=surfs;
      clear surfs
    end
  end
  if length(scale)==1
    shadeval=0.75;
  else
    shadeval=scale(2);
    scale=scale(1);
  end
  
  % Information about body and it's location (for axes and stuff)

  ld=numnodes*6;
  max(displacements);
  disx=displacements(1:6:ld,1);
  disy=displacements(2:6:ld,1);
  disz=displacements(3:6:ld,1);
  rots=[displacements(4:6:ld);displacements(5:6:ld);displacements(6:6:ld)];
  dismatrix=[disx disy disz];%,pause
  %dismatrix=dismatrix/max(abs(rots))*maxrotplot;
  maxdismag=max(max(abs(dismatrix)));%pause
  dismatrix=dismatrix/maxdismag*maxdisplot*maxnodespread;
  max(max(dismatrix));
%  if maxdismag>maxdisplot*maxnodespread
%    dismatrix=dismatrix/maxdismag*maxdisplot*maxnodespread;
%  end
  if maxdismag==0
    maxdismag=1;
  end
%  dismatrix=dismatrix*100
%  scale

  dnodes=nodes+scale*real(dismatrix);
  if exist('surfs')
    numsurfs=size(surfs,1);
  else 
    numsurfs=0;
  end
  %for nodenum=1:numnodes
  defnodeshandle=plot3(dnodes(:,1),dnodes(:,2),dnodes(:,3),'o');
  set(defnodeshandle,'markersize',4)
  cameratoolbar('show')
  hold on
  undefnodeshandle=plot3(nodes(:,1),nodes(:,2),nodes(:,3),'o');
  set(undefnodeshandle,'markersize',4,'color',[1 1 1]*shadeval)
  %numsurfs
  for surfnum=1:numsurfs
    coords=nodes(surfs(surfnum,1:4),:);
    patchcolor=surfs(surfnum,5:7);
    h=patch(coords(:,1),coords(:,2),coords(:,3),patchcolor+.7*([1 1 ...
		    1]-patchcolor));%patchcolor)
    alpha(h,.4)
  end
  for surfnum=1:numsurfs
    coords=real(dnodes(surfs(surfnum,1:4),:));
    patchcolor=surfs(surfnum,5:7);
    h=patch(coords(:,1),coords(:,2),coords(:,3),patchcolor* ...
	    shadeval);
    alpha(h,.8)
  end
  for linenum=1:numlines
    curlineu=real([nodes(lines(linenum,1),:);nodes(lines(linenum,2),:)]);
    undeflinehand=plot3(curlineu(:,1),curlineu(:,2),curlineu(:,3),'y');
    set(undeflinehand,'LineWidth',.5,'color',[1 1 1]*shadeval)
    curline=real([dnodes(lines(linenum,1),:);dnodes(lines(linenum,2),:)]);
    deflinehand=plot3(curline(:,1),curline(:,2),curline(:,3));
    set(deflinehand,'LineWidth',2)
  end
  %axis square
  hold off
  grid on
end
xlabel('x','buttondownf','edtext')
ylabel('y','buttondownf','edtext')
zlabel('z','buttondownf','edtext')
%maxx=max(disx);minx=min(disx);  
%maxy=max(disy);minx=min(disy);  
%maxz=max(disz);minx=min(disz);
%allnodes=[dnodes; nodes];
%maxspread=max([max(allnodes(:,1))-min(allnodes(:,1)) max(allnodes(:,2))- ...
%	       min(allnodes(:,2)) max(allnodes(:,3))-min(dnodes(:,3)) ]);
maxspread=maxnodespread*(1+.1+maxdisplot);
%maxspread=max([(maxsizex-minsizex)+abs(maxx),(maxsizey-minsizey)+ ...
%               abs(maxy),(maxsizez-minsizez)+abs(maxz)]) 
axis([-1 1 -1 1 -1 1]*maxspread/2*1.1+[avex avex avey avey avez ...
		    avez])
%axis
%T
if exist('T')
  AZ=T(1);EL=T(2);
  view(AZ,EL)
end
drawnow
