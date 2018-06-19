clc; close all; clear all;

%% To find Transfer Function for all five cases for mdofcf %
load frf1.mat
h1=frf1;
load frf2.mat
h2=frf2;
load frf3.mat
h3=frf3;
load frf4.mat
h4=frf4;
load frf5.mat
h5=frf5;
load frf6.mat
h6=frf6;
load frf7.mat
h7=frf7;
load frf8.mat
h8=frf8;
load frf9.mat
h9=frf9;
load frf10.mat
h10=frf10;
load frf11.mat
h11=frf11;
load frf12.mat
h12=frf12;
load frf13.mat
h13=frf13;
load frf14.mat
h14=frf14;
load frf15.mat
h15=frf15;
load frf16.mat
h16=frf16;
load frf17.mat
h17=frf17;
load frf18.mat
h18=frf18;
 
tf=[h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18];

%% Computing damping, natural frequency, mode shape vector  use mdofcf.m in vibration tool box

display('%%%%%%%% Computing damping, natural frequency, mode shape vector  use mdofcf %%%%%%%')
lamda=zeros(18); z=zeros(18); s=zeros(18);

f1 = 13.52; f2 = 51.61; f3 = 124.84; f4 = 167.18; f5 = 227.76;  % select the peak values from FRF plots for each case %
load freqencies.mat
%load Lab4_Case1;
figure
H=20*log10(abs(frf3)); % in db
plot(freqencies,H,'r'); hold on;
[z1,nf1,u1]=vtb7_4(freqencies,tf,f1-5,f1+5)   % use -5,+5 lower and upper bounds at the peak value
z(1,1)=z1;
lamda(1,1)=(2*pi*nf1)^2;
s(:,1)=real(u1)

% load Lab4_Case2;
H=20*log10(abs(frf1)); % in db
plot(Freq_domain,H, 'c'); hold on;
[z2,nf2,u2]=mdofcf(Freq_domain,tf,f2-5,f2+5);
z(2,2)=z2;
lamda(2,2)=(2*pi*nf2)^2;
s(:,2)=real(u2)

load Lab4_Case3;
H=20*log10(abs(Hf_chan_2)); % in db
plot(Freq_domain,H, 'b'); hold on;
[z3,nf3,u3]=mdofcf(Freq_domain,tf,f3-5,f3+5);
z(3,3)=z3;
lamda(3,3)=(2*pi*nf3)^2;
s(:,3)=real(u3)

load Lab4_Case4;
H=20*log10(abs(Hf_chan_2)); % in db
plot(Freq_domain,H, 'm'); hold on;
[z4,nf4,u4]=mdofcf(Freq_domain,tf,f4-5,f4+5);
z(4,4)=z4;
lamda(4,4)=(2*pi*nf4)^2;
s(:,4)=real(u4)

load Lab4_Case5;
H=20*log10(abs(Hf_chan_2)); % in db
plot(Freq_domain,H, 'g'); grid on; xlabel('Frequency - Hz'); ylabel('FRF - db'); hold on;
title('FRF vs Frequency')
legend('Case 1','Case 2', 'Case 3','Case 4','Case 5')
[z5,nf5,u5]=mdofcf(Freq_domain,tf,f5-4,f5+4);
z(5,5)=z5
lamda(5,5)=(2*pi*nf5)^2
s(:,5)=real(u5)

Dampingratio_mdofcf = [z1, z2, z3, z4, z5]            % Damping ratio's using mdofcf
Naturalfreq_mdofcf = [nf1, nf2, nf3, nf4, nf5]        % Natural Frequency using mdofcf
Modeshape_mdofcf = [u1, u2, u3, u4, u5]               % Mode Shape vector using mdofcf

m=(s'.*eye(5)).*s;              % Mass Matrix
c=(s'.*(2*z*sqrt(lamda))).*s;   % Damping Matrix
k=(s'.*lamda).*s;               % Stiffness Matrix

MM=vpa(m,8)  % use Variable Precision Arithmatic upto 5~8 %
CM=vpa(c,8)  % MM, DM, SM matices in the command window %
KM=vpa(k,8)  
