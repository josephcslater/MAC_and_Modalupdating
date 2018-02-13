%Some commonly used units that are not kms. Multiplying by these
%factors will convert to appropriate kms units. They may be
%compounded. 10*psi will give the converted value of 6.894757e4
%Pa. This could have aldo been done the hard way with
%10*lbf/in^2. See the actual code for currently implemented units.

%All values obtained from Mark's Handbook
in=2.54e-2;
sqin=in^2;
mm=1e-3;
cm=10*mm;

psi=6.894757e3;
ft=3.048e-1;
yard=3*ft;
slug=14.5939;
lbm=0.4535924;
lbf=4.448222;
g=9.80665;%  Standard gravity is defined to be that at sea level
          %  and latitude of 45 degrees. 
gm=1e-3;
angstrom=1e-10;
pm=1e-12;
nm=1e-9;
mum=1e-6;
deg=pi/180;
km=1000;

% E G rho alpha
aluminum6061=[ 10.6e6*psi 10.6e6*psi/(2*(1+.33)) 2770 (13.5e-6)*9/5];
aluminum=aluminum6061;
steelsheet=[  28e6*psi 28e6*psi/(2*(1+.29))  8300];
steelhot=[ 29.5e6*psi 29.5e6*psi/(2*(1+.29)) 7850];
steel=steelhot;