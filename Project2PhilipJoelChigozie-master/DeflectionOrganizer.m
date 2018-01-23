function DeflectionOrganizer( filename1,filename2,filename3 )
% DeflectionOrganizer takes input node locations and processed deflections
% and organizes them into vectors for plotting
% filename1 and filename2 can be in any order but must contain the
% variables nodes and Xg
% filename3 is the name of the file the vector data is written to
load(filename1);
load(filename2);
X=nodes(:,1);
Y=nodes(:,2);
Z=nodes(:,3);
U=Xg(:,1);
V=Xg(:,2);
W=Xg(:,3);
save(filename3,'X','Y','Z','U','V','W');

end

