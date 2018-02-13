function Eb = YoungsModulus3D(E, nu, G)
% Youngs Modulus Matrix in 3D based on Cook's Textbook 
% Page 79, Equation 3.1-5
c = E/(1+nu)/(1-2*nu);
E11 = (1-nu)*c;
E22 = E11;
E33 = E11;

E44 = G;
E55 = G;
E66 = G;

E12 = nu*c;
E21 = E12;
E13 = E12;
E31 = E12;
E23 = E12;
E32 = E12;

Eb = [E11, E12, E13, 0  , 0  , 0
      E21, E22, E23, 0  , 0  , 0
      E31, E32, E33, 0  , 0  , 0
      0  , 0  , 0  , E44, 0  , 0
      0  , 0  , 0  , 0  , E55, 0
      0  , 0  , 0  , 0  , 0  , E66];
