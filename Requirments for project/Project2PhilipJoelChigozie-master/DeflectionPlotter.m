function DeflectionPlotter( filename )
% DeflectionPlotter plots the points and vectors defined in filename on one
% plot
% filename is a .mat file containing variables X, Y, Z, U, V, and W
load(filename);
figure
hold on
scatter3(X,Y,Z,'filled');
quiver3(X,Y,Z,U,V,W);
hold off

end