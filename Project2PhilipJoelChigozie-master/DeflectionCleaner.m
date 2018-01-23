function DeflectionCleaner( filename,constant,exponent,sigfigs )
% DeflectionCleaner cleans output WFEM deflections with a lowpass
% magntitude filter and significant figure reduction. It returns a 3 x node
% count array with nodal deflections in it and overwrites the input file
% with it.
%
% filename is the name of the file containing the global deflection matrix
% 'Xg'. This file is a .mat file.
% constant is a double value used to lowpass filter 'Xg'
% exponent is the order of magnitude of the lowpass filter
% (constant*10^-exponent)
% if exponent is 0 the lowpass filter uses (constant)
% sigfigs is the number of significant figures to use
%
% Format example 1: DeflectionCleaner('filetoopen.mat',1.5,10,5)
% loads filetoopen.mat with variable Xg
% lowpass filters with cutoff of 1.5*10^-10
% rounds to 5 significant figures
%
% Format example 2: DeflectionCleaner('filewithXginit.mat',(1/1000),0,3)
% loads filewithXginit.mat with variable Xg
% lowpass filters with cutoff of 1/1000
% rounds to 3 significant figures

load(filename);
% Check to see if Xg is sparse
if issparse(Xg)
    Xg=nonzeros(Xg);
end
% set pseudo-zero values to zero
if exponent~=0
    for i=1:size(Xg,1)
        for j=1:size(Xg,2)
            if abs(Xg(i,j))<constant*10^-exponent
                Xg(i,j)=0;
            end
        end
    end
else
    for i=1:size(Xg,1)
        for j=1:size(Xg,2)
            if abs(Xg(i,j))<constant
                Xg(i,j)=0;
            end
        end
    end
end

% round other values to more reasonable precision
Xg=round(Xg,sigfigs,'significant');
% check if 1D array of values
if size(Xg,2)==1
    j=0;
    k=1;
    for i=1:size(Xg,1)
        j=j+1;
        Xgnew(k,j)=Xg(i);
        if j==3
            j=0;
            k=k+1;
        end
    end
    Xg=Xgnew;
end
save(filename,'Xg');
end
