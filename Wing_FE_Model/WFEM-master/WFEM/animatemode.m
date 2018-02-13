function animatemode(nodes,lines,surfs,X,T)

displacements=X;  
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
curtitle=get(get(gca,'title'),'string');
numsteps=20;
ld=length(displacements);
disx=displacements(1:6:ld,1);
disy=displacements(2:6:ld,1);
disz=displacements(3:6:ld,1);
dismatrix=[disx disy disz];
X=X/max(max(abs(dismatrix)))*0.1*maxnodespread;


%plot(X),pause
for phaseinc=1:numsteps
  cosphase=cos((phaseinc-1)/numsteps*2*pi);
  %pause
  plotgeom('deformed',nodes,lines,surfs,X*1.00,[ cosphase .75 ],T)
  %disp('hello'),pause
  title(curtitle)
  drawnow
  figure(gcf)
    %plotcommand
  M(phaseinc)=getframe;
% need avifile command before this
%  MOV=addframe(MOV,M(phaseinc));
end
movie(M,5);
answer=lower(input('Do you want to play the animation again? ', ...
		   's'));
while length(findstr(answer,'y'))>0
  movie(M,5);
  answer=lower(input('Do you want to play the animation again? ', ...
		   's'));
end

  
answer=lower(input('Do you want to save the movie as an AVI file? ', ...
                   's'));
if length(findstr(answer,'y'))>0
  [avifilename,pathname]=uiputfile('*.avi','Save AVI movie as:');
  while length(avifilename)==1
    [avifilename,pathname]=uiputfile('*.avi',['Save AVI movie as (put '...
                                     'name in lower box);']);
  end
  if isstr(avifilename) & isstr(pathname)
    curdir=pwd;
    cd(pathname)
    avihandle=avifile(avifilename);
    for i=1:numsteps
      avihandle=addframe(avihandle,M(i));
    end
    avihandle=close(avihandle);  
    cd(curdir)
  end
end

%movie2avi(MOV,filename)
% Need to be able to repeat movie.