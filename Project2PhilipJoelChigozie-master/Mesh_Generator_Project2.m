clear all
close all
clc

Name = 'SquareBeam';
filename=sprintf('%s.txt',Name); %writes the input file name
fid=fopen(filename, 'w'); %Opens the input file for writing

%% 3-D brick mesh.

xLength = 74; xDelta = 18.5; xArray = -xLength/2:xDelta:xLength/2;
yLength = 74; yDelta = 18.5; yArray = -yLength/2:yDelta:yLength/2;
zLength = 190; zDelta = 38; zArray = 0:zDelta:zLength;

xLengthTop=74;
yLengthTop=74;
numz=size(zArray,2);
draftx=linspace(1,(xLengthTop/xLength),numz);
drafty=linspace(1,(yLengthTop/yLength),numz);

[X,Y] = meshgrid(xArray,yArray);
[row,col] = size(X);

numberOfNodes = length(zArray)*row*col;
nodes = zeros(numberOfNodes,4);

NumberOfNodesPerLevel = numberOfNodes/length(zArray);
NumberOfElements = (length(xArray)-1)*(length(yArray)-1)*(length(zArray)-1);

index = 1;
for k = 1:length(zArray)
    for i = 1:row
        for j = 1:col
            nodes(index,:) = [index , X(i,j)*draftx(k), Y(i,j)*drafty(k), zArray(k)];
            index = index +1;
        end
    end
end

figure, hold on
scatter3(nodes(:,2),nodes(:,3),nodes(:,4),150)
[rn,cn]=size(nodes);

for i = 1:numberOfNodes;
   text( nodes(i,2),nodes(i,3),nodes(i,4),int2str(nodes(i,1)))
end

nx = length(xArray); % # nodes in the X
ny = length(yArray); % # nodes in the Y

elementConnectivity = [(1:NumberOfNodesPerLevel)', (2:NumberOfNodesPerLevel+1)', (nx+2:NumberOfNodesPerLevel+nx+1)', (nx+1:NumberOfNodesPerLevel+nx)'];

ghostElementsLeft = (nx):(nx):(nx)*(ny);
ghostElemetsTop = (ny-1)*nx:nx*ny;
removeElemetns = union(ghostElementsLeft, ghostElemetsTop);

elementConnectivity(removeElemetns, :) = [];
elementConnectivity

connections = [];

for i = 1:length(zArray)-1
    plusNumberOfNodesPerLevel = elementConnectivity+NumberOfNodesPerLevel*i;
    connections = [connections; [elementConnectivity+NumberOfNodesPerLevel*(i-1), plusNumberOfNodesPerLevel] ];
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~Printing to a file~~~~~~~~~~~~~~~~~~~~~~~~%

fprintf(fid,'%% variables\n');
fprintf(fid,'%% All of these actions are not the most efficient for this problem.\n');

fprintf(fid,'\n');

fprintf(fid,'element properties\n');
fprintf(fid,'%% Beam format\n');
fprintf(fid,'%% E G rho\n');
fprintf(fid,'steel(1:3)\n');

fprintf(fid,'\n');

fprintf(fid,'BrickNoCorrections elements\n');
fprintf(fid,'%%node1 node2 node3 node4 node5 node6 node7 node8 propertynumber points\n');

[rcon,ccon]=size(connections);

for i = 1:NumberOfElements
    
     fprintf(fid,'%d %d %d %d %d %d %d %d 1\n',connections(i,:,:,:,:,:,:,:));
 
end

fprintf(fid,'\n');

fprintf(fid,'nodes\n');
fprintf(fid,'%% I can include comment lines\n');
fprintf(fid,'%%node num, x y z, Node number isnt ever stored in nodes matrix\n');

for i = 1:rn
  
    fprintf(fid, '%d %d %d %d\n',nodes(i,:,:,:));
  
end

fprintf(fid,'\n');

fprintf(fid,'points\n');
fprintf(fid,'1 1 1 1\n');

fprintf(fid,'\n');

fprintf(fid,'fix clamp\n');
for i=1:NumberOfNodesPerLevel
    fprintf(fid,'%i\n',i);
end

fprintf(fid,'\n');

fprintf(fid,'staticanalysis\n');
fprintf(fid,'plotdeformed\n');
fprintf(fid,'X\n');

fclose(fid);