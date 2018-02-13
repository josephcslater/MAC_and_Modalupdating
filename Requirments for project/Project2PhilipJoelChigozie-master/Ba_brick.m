function [ Ba ] = Ba_brick( J,Epsilon,Eta,Zeta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Ba1=J\[-2*Epsilon; 0; 0];
Ba2=J\[0; -2*Eta; 0];
Ba3=J\[0; 0; -2*Zeta];
Bi1=[Ba1(1) 0      0;
     0      Ba1(2) 0;
     0      0      Ba1(3);
     Ba1(2) Ba1(1) 0;
     0      Ba1(3) Ba1(2);
     Ba1(3) 0      Ba1(1)];
 Bi2=[Ba2(1) 0      0;
     0      Ba2(2) 0;
     0      0      Ba2(3);
     Ba2(2) Ba2(1) 0;
     0      Ba2(3) Ba2(2);
     Ba2(3) 0      Ba2(1)];
 Bi3=[Ba3(1) 0      0;
     0      Ba3(2) 0;
     0      0      Ba3(3);
     Ba3(2) Ba3(1) 0;
     0      Ba3(3) Ba3(2);
     Ba3(3) 0      Ba3(1)];
 Ba=[Bi1 Bi2 Bi3];

end

