Globals
---------

=============     ===========
Variable name     Description
=============     ===========
M	                Global Mass matrix
K                 Global Stiffness Matrix
Bso               Second order form input matrix
Bthermal          Second order form thermal input matrix.
=============     ===========



Locals
------------

=============     ===========
Variable name     Description
=============     ===========
Mr                Reduced Mass matrix, T	T'MT gives Mr
Kr                Reduced Stiffness matrix, T'KT gives Kr
x                 displacements, x=T*xr
x                 reduced displacements
slavedofs         Slave DOFs
=============     ===========
