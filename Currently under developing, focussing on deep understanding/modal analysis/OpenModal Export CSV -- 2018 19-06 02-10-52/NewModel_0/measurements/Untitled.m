

realpart= csvread('measurement_17.csv',1,1,[1 1 4097 1]);
     imagpart= csvread('measurement_17.csv',0,1,[2 2 4097 2]);
     imagpart=imagpart(:,2)
     frf18 = realpart + (imagpart)*i

save frf18.mat
%save freqenci

plot(freqencies,20*log10(frf18))
title('frf18')