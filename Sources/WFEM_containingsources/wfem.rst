WFEM
######
This manual isn't quite translating properly. Please view WFEM.pdf for now. 


Wright Finite Element Method (WFEM) is designed to be an easily
extensible research-oriented structural finite element code. The design
concept is that all “parts” of the model are actions by external codes.
The central finite element code knows nothing about any particular
element or boundary condition type. It depends on a completely
self-capable code for implementing the “method” of the element or
boundary condition. The advantage of this approach is that the code is
able to handle the introduction of new elements without any modification
of existing parts of the package. Further, additional capabilities can
be added to elements without disruption to the main body of the code.

Anatomy of the body of the code
===============================

The main body of FEM is currently under construction, and likely will be
indefinitely. In terms of performing analysis, this is not all that bad,
since the real work of the code is handled internally to the element
codes. The ability to find mode shapes and natural frequencies, plot
deflections, and find stresses are separate aspects from the main body,
and are performed by the element codes.

Running WFEM
------------

To run WFEM, set your M\ atlab path to the WFEM directory and type .
WFEM will prompt you for an input file name. Alternatively, type for
other means by which to start a run of WFEM.

Linear Dynamic Analysis
-----------------------

Finding mode shapes and natural frequencies of a linear model can be
performed using which uses Cholesky decomposition of the mass matrix and
assumes nothing about the stiffness matrix, :math:`K`, in terms of
symmetry. However, no asymmetric capabilities have been added yet (the
FEM code doesn’t have the ability to make an asymmetric element).

Guyan Reduction
---------------

Guyan reduction of the model can be easily performed using . In the
automatic mode, it will remove DOFs with a mass/stiffness below a
predefined threshold.

Plotting
--------

The ability to be plotted is one that is defined internally to an
element. The plotting routine knows nothing about elements, only about
nodes, their locations, and how they are connected (via lines,
surfaces...). Lines and patches are defined by what nodes they connect,
not by geometry.

To draw lines within elements, see the element. Look for the lines
variable. Appending the pair of node numbers to be connected by line to
this array ensures that it will be drawn by WFEM. Explanation for how to
draw a surface is also given in the comments.

To draw surfaces, see the generate section of the element. The color of
a quadrilateral patch is defined by . Adding a patch to be drawn by
listing its node numbers (clockwise or counterclockwise) followed by its
color is all that’s needed to assure that it is drawn in static and
dynamic situations.

Plotting of modes is easily accomplished internally by WFEM by adding
mode shapes to the undeformed nodal location, and plotting all
predefined nodes, lines, and patches. Element authors need do nothing to
enable plotting of mode shapes/deformed static shapes.

Plotting Stresses/Strains (Not completed)
-----------------------------------------

Plotting stresses and strains is accomplishes by using the and matrices.
These are actually stresses and strain levels corresponding to the
appropriate lines and patches and are placed in those matrices at the
correct locations by the appropriate element code.

Entering data
=============

For the sake of ease of use, it is best to define nodes and element
property types before elements are defined. An example input code is
shown in listing [listing:example.txt] on page . The first thing that
one notices about the file is that the ’%’ can be used to add comments
throughout the file, at any location, even in the middle of data. The
’%’ must reside in the *first column* for the comment to be recognized
as a comment character. That means *you cannot put a comment at the end
of a line of data.*

Properties are defined after the ‘element properties’ line. Some
flexibility is allowed for typographical errors, but this flexibility
shouldn’t be considered a poetic license. Each type of element can have
its own method for defining element properties. Elements can also have
multiple methods of data input. For example, it would be quite
reasonable to allow for a beam element that has constant properties
along its length. As a result, only 7 entries,
:math:`\rho, A, I_{yy}, I_{zz}, J, E,` and :math:`G` would be required.
For example, see the element documentation in section [el:beam3]
regarding allowable material input formats and how they are interpreted.
The WFEM code makes no attempt to understand these values at all. It
simply loads them and stores them where the element codes can easily
find them. Element property indices are identified by their row number.
In the example shown (See Listing [listing:example.txt]), there are two
element properties given. All elements defined later must refer to one
of these two definitions unless the element has no properties (we’re
staying very flexible here.).

Nodes are defined subsequent to the command (see section
[command:nodes]). Unlike elements, the first value on each line is the
node number. This example shows some of the flexibility of the reader in
that it allows nodes to be defined before or after the elements. In
fact, you can alternately define nodes and elements. This is convenient
if you have an existing model and would like to simply tack on an
additional substructure to the end.

Elements are defined using ‘*elementname* ’ (in this case ). There must
exist a function *elementname*. Rows subsequent to such a command are
defined to be of that type. The definition of elements is concluded with
a blank line. Keep in mind that the main body of the finite element code
knows *nothing* about any element, so adding an element to the code
means writing a function file to do all of the work of using the
element.

Units
-----

Like any code, as long as consistent units are used, the code doesn’t
care what units you are using. The presumption of the code, however, is
that values are entered in mks units (kilogram, meter, seconds). To
assist with this, some built in constants are available for doing
conversions for the user. For example, to enter 1 inch using presuming
mks units, one would instead type *1\*in*. Conversions from mks to th
inch-pound system are not available directly, but can always be achieved
by inverting the conversions. Currently defined units are *in, ft, yard,
lbf, lbm, slug, psi, deg, gm, pm, angstrom, nm, :math:`\mu`\ m (mum),
mm, cm* and *km*.

Defined Constants
-----------------

The gravitational constant *g* is also defined to be 9.80665 (standard
gravity is defined to be at sea level and latitude
45\ :sup:`:math:`\circ``.

Material propeties, E, G, :math:`\rho`, are available instead of typing
them out as the variables *aluminum6061* (also *aluminum*),
*steelsheet,* and *steelhot* (also *steel*). Others can easily be added
to .

Variables
---------

| Variables can also be used to define geometries. Variables are
  generated by definitions (e.g. ) following the variable command ( on a
  line by itself). Variable can be defined using any lowercase letter
  with the exceptions of :math:`i` and :math:`j` which are reserved for
  :math:`i=j=\sqrt{-1}`. Variables can also be lowercase letters
  followed by a number. Other variable names may be used with care. The
  guidelines given here are provided to simplify selection of names. As
  long as the variable name doesn’t conflict with an internally used
  variable name, it is acceptable to use. A warning will be generated if
  an illegal choice is used.
| Example:

::

      command
      s=1

Parameters
----------

The third, optional, argument to is a vector of parameters that are used
in the problem definition of the input file. If this vector is
available, the input file may use variables *parameter(i)* where *i* has
a numerical value between :math:`1` and the total length of the input
parameter vector. The advantage of parameters is that multiple similar
cases can be run without having to repeatedly edit the input file. For
instance,

::

      for i=1:10
      wfem('inputfilename.txt',[],i)
      end

would run the case defined in to be run 10 times, with the quantity
referred to as *parameter(1)* in the input file varying between 1 and
10.

Defining Nodes
--------------

Nodes are defined subsequent to the command. the format is

::

    nodes
    nodenum xlocation ylocation zlocation
    %This is a comment

Definition of nodes ends with a blank line. Note that the node numbers
are listed in the format for convenience only. The actually assigned
node number is incremented by one for each new line of data. That is,
node numbers monotonically increase by one.

Meshing
=======

| Meshing in this code is different than meshing in a traditional finite
  element code. In traditional finite elements, a geometry of an object
  is defined, the internal region of the body is divided by algorithm
  into elements, then nodes are automatically generated. In WFEM,
  meshing is the repeating of a substructure multiple times in a one
  dimensional arrangement. To do this, the command is entered.
  Everything between this command and the command is considered to be
  the definition of what one bay looks like. This should include *only*
  elements. Elements must be *bay capable*\ (See section
  [sec:elements]), as the bay meshing requires some knowledge to be
  obtained regarding the elements. Including non-bay capable elements
  will lead to erroneous results. No error checking for this exist. The
  command is formated as:
| *l* *M* *N*.
| Here the nodes listed in *M* are the innermost nodes of the unit bay,
  and the nodes listed in *N* are the outermost nodes, or free nodes.
  Another way to think of this is that the nodes *N* are currently the
  furthest extending nodes in what will become a series of bays. The
  nodes *M* are the starting nodes of this unit. For the next bay, the
  starting nodes will be *N*.For example, the code segment:

::

    bay element
    beam3 element
    1 2 1

    repeat bay 3 times. Attach 1 to 2.

will create will create a four element model of a beam by creating the
first element (connecting nodes 1 and 2) then creating 3 more elements.

Here *l* defines how many times the bay is to be repeated. *M* is a list
of nodes defining the attachment nodes of the bay that are attached to
the existing defined structure. These can be thought of as the starting
nodes of the bay. These node numbers will correspond to nodes *N* of the
next bay, where *N* is the list of the end, or exposed, nodes of the bay
that subsequent bays will be attached to.

Boundary Conditions and Constraints
===================================

Boundary conditions and constraints are applied using the and commands.
They are coupled with the appropriate entity as *clamp* or *pin*. The
former command rigidly connects two nodes in all six degrees of freedom,
and the second results in a pin connection at subsequently listed nodes.
The commands are followed by lines describing individual boundary
conditions/constraints, the format being specific to the type of
constraint. The following table describes all available attachment
types.

+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| Attachment name   | Description                                 | Input Format                              | NASTRAN Support   |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| clamp             | all degrees of freedom                      | node1 (node2), nodes must be coincident   | Y                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| pin               | Thrust bearing hinge                        | node1 (node2) udv                         | N                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| roller            | Roller in track                             | node1 (node2) udv hinge, udv motion       | N                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| ball              | Ball joint                                  | node1 (node2)                             | N                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| rod               | rod with pins at each end                   | node1 node2                               | N                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| rbeam             | rigid (clamped) beam                        | node1 node2                               | Y                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| surfaceball       | 1 translation restricted                    | node1 (node2) udv translation             | N                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| surface           | 1 trans + 3 rots restricted                 | node1 (node2) udv translation             | N                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+
| rigidbody         | calculate inertia tensor relative to node   | node1 0                                   | N                 |
+-------------------+---------------------------------------------+-------------------------------------------+-------------------+

#. Unit direction vectors must be orthogonal.

In cases of constraints, the appropriate DOFs of the second node are
reduced (slaved) to those of node 1.

The final attachment name is a constraint condition where rigid massless
beams are presumed to connect all nodes, reducing the model to the
single node *node1*. All other constraints will be ignored under these
circumstances. The dummy value :math:`0` is inserted for compatibility
with other constrain condition formats. The action is necessary to
obtain the rigid body parameters. See the file for example usage.

See example.txt (Listing [listing:example.txt]) for an example
application of a boundary condition.

Applied Static Loads
====================

Application of a static load is performed using the command. Subsequent
lines list the node number, the degree of freedom number, and the amount
of the load applied to that node at that degree of freedom. Loads must
be in consistent units. See Listing [listing:example.txt] for an example
of load application.

Actions
=======

Just as element routines act on the modal, more global actions are
performed by action routines. Many actions are performed by default when
a model is built. These include , , , and , their function being clear
by their names. Additional detail regarding these action can be gleaned
from the help provided in the code (for example ). Additional actions
currently available follow. Actions are intelligent enough to recognize
when they cannot be performed without another preceding analysis. In the
case that another action must take place first, an error is produced
stating explicitly what action was needed.

plotundeformed
--------------

The action plots the as-defined structure and geometry, a useful tool in
finding date entry errors.

findinitialstrain
-----------------

The action will apply the initial strain inducing loads and return the
initial global deflections in the variable *X*. The initially deformed
structure is automatically drawn if graphical output is available.

staticanalysis
--------------

The action will apply the initial strain inducing loads and the
prescribed loads and return the global deflections in the variable *X*.
The deformed structure is automatically drawn if graphical output is
available.

plotdeformed
------------

The action will plot the deformed structure under the previously solved
for deflections *X*.

modalanalysis
-------------

The action obtains the natural frequencies (in Hz) and mode shapes and
stored them in the variables *fs* and *fms*. These results are also
automatically saved to the restart file.

modalreview
-----------

The action is a poor-man script that runs through the results of modal
analysis. Not very sophisticated, but it gets the job done.

reducedofs
----------

The action will reduce the degrees of freedom in the system to remove
DOFS associated with boundary conditions and constraints. Matrices
returned are *Mr*, *Kr* and *Tr* where

.. math:: M_r=T_r^T M T_r

 and

.. math:: K_r=T_r^T K T_r

Find Inertia
------------

[totalmass,INERTIATENSOR,CG]=FINDINERTIA solves for the inertia tensor
relative to the center of gravity, as well as the center of gravity.

planefit
--------

The action fit the best possible plane to the nodes containing surface
elements. For best results, this plane should b closer to perpendicular
to the :math:`z` axis than either of the others.

end
---

The action will cause wfem to be exited, leaving the user back are the
M\ atlab prompt.

M commands
----------

A great deal of flexibility is obtained by recognizing that any M\ atlab
command can be executed as an action. This allows manipulation of
graphics, intermediate computations, and other miscellaneous actions not
listed here to be performed. A warning will be given that the command
are not explicitly defined, and that they are being passed directly to
the M\ atlab interpreter.

Elements
========

Element modes
-------------

Elements are implemented through modal function calls to m-files of the
same name as the element. The mode is essentially “what action should be
taken”. Additional information is passed to the element as needed,
appropriate to the mode. Some modes must exist for WFEM to be able to do
anything with the element. Some are recommended, which is admittedly
subjective. Recommended means that the elements can be read into the
data structure, and I think those modes are not necessary for elementary
types of analysis (mode shapes, static deflection, linear simulation).
Optional modes are those such as dynamic loads due to rotating
coordinate systems, nonlinear stiffness terms, and stress and strain
calculation capability.

Required element modes
----------------------

takes the element data and puts all available element information into
the next available element entry in the element data structure. If
element information is missing, generate may optionally generate it
(i.e. internal nodes in are handled this way).

takes element properties and nodal locations, generates the local
coordinate finite element matrices, rotates them to global coordinates,
and assembles them into the stiffness and mass matrices. It also
generated the drawing properties for the element (points, lines, or
surfaces.). These are used for easy and faster drawing of the deformed
and undeformed structure.

is also used to generate entries in the control matrix, as well as
nonlinear force flags in the nonlinear force vector. These allow
generation of linear matrix models for control law design (in the first
case) as well as nonlinear time simulations.

Highly recommended modes
------------------------

Optional modes
--------------

numofnodes
~~~~~~~~~~

| The mode is executed as
| *elementname*\ *eldata*
| where *eldata* is the row vector defining the element. It returns as
  an output the number of the defining parameters that are node numbers,
  starting from the left. This is required for , and thus is optional if
  the element isn’t to be used for meshing.

Element tools available
-----------------------

Gauss Legendre Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~

The file makes Gauss Legendre integration in 1-3 dimensions and up to
10th order simple with a single function call and a loop. Weights and
gauss points are obtained using *pg*\ *wg*\ *n* where *n* is a vector of
length 1-3, each entry determining how many integration points to use in
that direction. *pg* and *wg* are the resulting Gauss points and
weights. *pg* is in the form

.. math::

   pg=\begin{bmatrix}
   x_{1}&y_{1}&z_{1}\\
   x_{2}&y_{2}&z_{2}\\
   \vdots&\vdots&\vdots\\
   \end{bmatrix}

\ and *wg* is simply a list of weights corresponding to the points of
*pg*.

Element types available
-----------------------

The element is a simple constant-cross section two-noded rod (spring)
element with a consistent mass matrix. It has linear shape functions.
Cross section properties are assumed to be constant.

| The element may be defined by any of the following:

::

      rod3 elements
      node1 node2 node3 pointnumber materialnumber
      node1 node2 pointnumber materialnumber
      node1 node2 materialnumber

Only the first two nodes and the material number matter. The remainder
of the input options are available only for simplicity of switching from
elements to elements. See section [el:beam3] on elements for details of
the meaning of other values.

Material properties for the element are better referred to as rod
properties. They are entered using the material properties command and
can be entered as any of the following forms:

::

      element properties
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3
      E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2
      E G rho A J Izz Iyy
      E G rho A J Izz Iyy sx2 sy2 sz2 srx2 sry2 srz2 distype
       E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          sx2 sy2 sz2 srx2 sry2 srz2 distype
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 distype
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz3 mx3 my3 mz3 mrx3 mry3 mrz3
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...
          mx3 my3 mz3 mrx3 mry3 mrz3 sx3 sy3 sz3 srx3 sry3 srz3 distype
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...
          lx2 ly2 lz2 lrx2 lry2 lrz2
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...
          lx2 ly2 lz2 lrx2 lry2 lrz2 ...
          mx3 my3 mz3 mrx3 mry3 mrz3 sx3 sy3 sz3 srx3 sry3 srz3 ...
          lx3 ly3 lz3 lrx3 lry3 lrz3
      E rho A
      E rho A m
      E rho A m s
      E rho A m s distype
      E rho A m l distype
      E rho A m s l

Property formats other than the last 6 add no capabilities but are
available only for substituting elements for elements. , , and stand for
mean, standard deviation, or limit. A of :math:`0` means Gaussian
distribution. A of :math:`-1` means uniform distribution.

The element is a simple constant-cross section two-noded rod (spring)
element with a consistent mass matrix. It has linear shape functions.
Cross section properties are assumed to be constant. In can handle
initial strain definitions as well as thermal strains that change over
time.

| The element may be defined by any of the following:

::

      rod3 elements
      node1 node2 materialnumber

Material properties for the element are better referred to as rod
properties. They are entered using the material properties command and
can be entered as any of the following forms:

::

      element properties
      E rho A m s distype alpha alim thermalinput
      E rho A m l distype alpha alim thermalinput
      E rho A m s l       alpha alim thermalinput

The second 3 variables are for initial strains. , , and stand for mean,
standard deviation, or limit. A of :math:`0` means Gaussian
distribution. A of :math:`-1` means uniform distribution.

is the thermal coefficient of expansion. defines the maximum positive
deviation of , :math:`\left(\delta\alpha\right)_\text{max}`. allows
temperature inputs to be applied to multiple elements simultaneously.
:math:`\alpha` is allowed only to vary linearly at this time. The
thermal input matrix is , with columns corresponding to the numbers. The
matrix must be multiplied by the *deviation* from the nominal
temperature to obtain the psuedo-applied thermal load.

The element is a three-noded Euler-Bernoulli beam/rod/torsion element.
It has sixth order shape functions for the beam deformations, and
quadratic shape functions for torsion and extension deformations.
Properties vary quadratically (or lower order) along the length of the
element. The beam is asymmetric, allowing definition of the local
coordinates for each element (see below). Second moment of area
properties (:math:`I_{yy}` and :math:`I_{zz}`) must be defined in the
principle area coordinate frame. :math:`J` is the effective torsional
polar moment of area. It should be equal to :math:`I_{yy}+I_{zz}` for
circular cross sections, and less than the true polar moment of area for
non-circular cross sections. For torsional inertia of the element,
:math:`I_0` is calculated appropriately as :math:`I_0=I_{yy}+I_{zz}`.
Nodes are numbered 1-3-2, so that the extra node (3) is in the middle of
the beam. Node three must be placed precisely in the middle of the beam
if it is defined explicitly. Deviation will cause errors.

| The element may be defined by any of the following:

::

      beam3 elements
      node1 node2 node3 pointnumber materialnumber
      node1 node2 pointnumber materialnumber
      node1 node2 materialnumber

The point number is defined using the command (see listing
[listing:example.txt], page ) in the same fashion as the command
(Section [command:nodes]). The location of the point, along with nodes 1
and 2, defines the :math:`x-y` plane. The local :math:`x` axis is from
node 1 to node 2. The local :math:`y` axis is from node 1 to the point,
but only the component perpendicular to the :math:`x` axis. The
:math:`z` axis is then known via the right hand rule.

+----------------------------+-------------+-------------+-------------+-------------+
|                            | SS mode 1   | SS mode 2   | FF mode 2   | FF mode 3   |
+============================+=============+=============+=============+=============+
| Continuous theory          | 9.87        | 39.48       | 22.40       | 61.62       |
+----------------------------+-------------+-------------+-------------+-------------+
| Single 6th order Element   | 9.87        | 39.65       | 22.56       | 63.54       |
+----------------------------+-------------+-------------+-------------+-------------+
| 2 4th order Elements       | 9.91        | 43.82       | 22.42       | 70.17       |
+----------------------------+-------------+-------------+-------------+-------------+

Table: Quality of 6th order beam element for simply supported and
free-free frequency determination. Coefficient to
:math:`\sqrt{EI/\rho A l^{2}}`.

+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Kind                           | :math:`\omega=`                                | Boundary Cond.   | Mode #   | % Error           |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Rod (Extension)                | :math:`\frac{\pi}{l}\sqrt{\frac{E}{\rho}}`     | FF               | 2        | 10.27             |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Rod (Extension)                | :math:`\frac{\pi}{2l}\sqrt{\frac{E}{\rho}}`    | CF               | 1        | 0.375             |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Rod (Extension)                | :math:`\frac{3\pi}{2l}\sqrt{\frac{E}{\rho}}`   | CF               | 2        | 21.4              |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Rod (Torsion)                  | :math:`\frac{\pi}{l}\sqrt{\frac{G}{\rho}}`     | FF               | 2        | 10.27             |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Rod (Torsion)                  | :math:`\frac{\pi}{2l}\sqrt{\frac{G}{\rho}}`    | CF               | 1        | 0.375             |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Rod (Torsion)                  | :math:`\frac{3\pi}{2l}\sqrt{\frac{G}{\rho}}`   | CF               | 2        | 21.4              |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Beam (Local :math:`y` plane)   | :math:`22.3733\alpha`                          | FF               | 2        | 0.7               |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Beam (Local :math:`y` plane)   | :math:`61.6728\alpha`                          | FF               | 3        | 3.1               |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Beam (Local :math:`y` plane)   | :math:`\pi^2\alpha`                            | SS               | 1        | :math:`<`\ 0.05   |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Beam (Local :math:`y` plane)   | :math:`4\pi^2\alpha`                           | SS               | 2        | 0.4               |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Beam (Local :math:`z` plane)   | :math:`22.3733\alpha`                          | FF               | 2        | 0.7               |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Beam (Local :math:`z` plane)   | :math:`61.6728\alpha`                          | FF               | 3        | 3.1               |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Beam (Local :math:`z` plane)   | :math:`\pi^2\alpha`                            | SS               | 1        | :math:`<`\ 0.05   |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+
| Beam (Local :math:`z` plane)   | :math:`4\pi^2\alpha`                           | SS               | 2        | 0.4               |
+--------------------------------+------------------------------------------------+------------------+----------+-------------------+

Table: Single finite element natural frequency estimate error relative
to continuous theory. Circular cross sections assumed, with
:math:`G=1.5911\times 10^{10}`, :math:`E= 4.1369\times 10^{10}`,
:math:`l=3.048`, :math:`\rho=1.6608\times 10^3`,
:math:`I_{yy}=I_{zz}=7.0612\times 10^{-7}`, and
:math:`A=2.4322\times 10^{-4}`. For beam frequencies,
:math:`\alpha=\sqrt{\frac{E I}{\rho A l^4}}`.

Material properties for the element are better referred to as beam
properties. They are entered using the material properties command and
can be entered as any of the following forms:

::

      element properties
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3
      E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2
      E G rho A J Izz Iyy
      E G rho A J Izz Iyy sx2 sy2 sz2 srx2 sry2 srz2 distype
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          sx2 sy2 sz2 srx2 sry2 srz2 distype
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 distype
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz3 mx3 my3 mz3 mrx3 mry3 mrz3
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...
          mx3 my3 mz3 mrx3 mry3 mrz3 sx3 sy3 sz3 srx3 sry3 srz3 distype
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...
          lx2 ly2 lz2 lrx2 lry2 lrz2
      E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
          mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...
          lx2 ly2 lz2 lrx2 lry2 lrz2 ...
          mx3 my3 mz3 mrx3 mry3 mrz3 sx3 sy3 sz3 srx3 sry3 srz3 ...
          lx3 ly3 lz3 lrx3 lry3 lrz3

If the second format is used, a linear interpolation of the properties
is presumed. If the third format is used, constant properties are
assumed. In subsequent lines, initial deflections can be prescribed,
deterministically, or stochastically (random). The character stands for
*mean* value, and stands for *rotation* value. If a stochastic form is
used, a distribution type must be prescribes. A of means normal
distribution with standard deviation of values as prescribe. A of means
uniform distribution with bounds relative to the mean set by the values.
A truncated Gaussian distribution is demonstrated in the last line.
There the values are limits relative to the mean (just as the values for
a uniform distribution). Note that the dots illustrate continuation of
the same line. Continuations of a line using the M\ atlab notation …can
only be used in defining element properties. Torsional rigidity,
:math:`J`, must be less than or equal to :math:`Iyy+Izz` at any given
cross section.

In addition, the variable can be defined to override the warnings that
occur if the beam is defined to have :math:`\frac{D}{l}<0.1`. This is
generally not advisable, as beams with a small value of
:math:`\frac{D}{l}` tend to behave much more like Timoshenko beams.

The element adds point masses and inertias. The inertia properties can
be entered in the following forms and are treated as a “material type”.
In order of increasing complexity, the acceptable formats are:

::

      element properties
      m I
      m Ixx Ixy Ixz Iyy Iyz Izz
      m Ixx Ixy Ixz Iyy Iyz Izz l1x ly1 lz1 l2x l2y l2z

where the local :math:`x` coordinate points in the direction of the
vector [l1x l1y l1z], and the local :math:`y` coordinate points in the
direction of the vector [l2x l2y l2z]. Elements are defined using

::

      inertia elements
      node materialnumber

A element is a method used to mark panels for obtaining geometric output
later. Element properties are simply the *total* mass of the panel:

::

      element properties
      m

Elements are defined using

::

      inertia elements
      node1 node2 node3 node4 materialnumber

A element provides outputs of displacement, velocity, or acceleration to
. Property values 1, 2, and 3 represent displacement, velocity, or
acceleration respectively. All 6 DOF are measured. If you want fewer,
ignore the extras, or dispose of the appropriate rows in the output
matrices.

A element is a method used to mark nodes on a surface of importance for
any given reason. The set of surface nodes provides the points included
in any accuracy calculations to determine surface quality. have no
material property, and are defined by only the node numbers.
