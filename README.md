
This is a study to explore Modal Analysis, Correlation and Modal Updating of FE models using Test data. 

The project is systematically explained in the [Python jupyter notebook](https://github.com/sainag2473/MAC_and_Modalupdating/blob/master/MAC_JupyterNotebook/MACandModelCorrection.ipynb)

## Purpose of bridging the gap between Testing and FE-Simulation:


**Product Development Process:** There are five phases in this process along the cost curve.

- Top two phases: **Concept and Detail Drawing**, in these phases, the CAD and FEA simulation play the major role. So, any change of product design is very inexpensive and it doesn't take much time for modifications.


- Next two phases: **Prototype and Production**: in these phases, it becomes little bit more expensive to change the aspects of the product design. So, it gets more expensive as we move forward.


- Final phase: **Field Failure (Testing)**, this is the phase where product recalls occur to make a major design change and it gets much more expensive.

To avoid this sort of serious field failure, it is important to combine the testing and the simulation so that validating not only the physical design but also whether the FE simulation is correct. So, the first preference is given to validation at the test condition to use the right modeling technique in simulation to look at the variations and develop the model a head of time to see certain phenomenon that would replicate the real world.

Repository contains details of the following:
FE models using 6 d.o.f/node Beam element and 8 noded Brick element were coded in MATLAB which had been validated with Ansys Apdl. 

Testing was done in a laboratory environment using an open source experimental modal analysis software called "OpenModal" which had been validated with existing Bobcat commercial software.



## References to Modules and Books Used:

[WFEM](https://github.com/josephcslater/WFEM) : by Dr. Joseph Slater

[OpenModal](https://github.com/openmodal/OpenModal) : created as an Opensource software written in 
Python for Vibration Testing with excellent capabilities. 

[vibrationtesting](https://github.com/Vibration-Testing/vibrationtesting) : by Dr. Joseph Slater

[Vibration Toolbox]( https://github.com/vibrationtoolbox/vibration_toolbox.git) : by Dr. Joseph Slater
  
Finite Element Model Updating in Structural Dynamics by M.I. Friswell and J.E. Mottershead

Concepts and Applications of Finite Element Analysis, Fourth edition by Robert D. Cook, David S. Malkus, Michael E. Plesha and Robert J. Witt

[Siemens Corp. MAC](https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Modal-Assurance-Criterion-MAC/ta-p/368008)

Some Pictures and Bulletin points were taken directly from Modal Analysis Conference conducted by Bruel and Kjaer, Chicago, May2018.
Mr Robert Trepanier is the presenter who has more than 30 years of experience, I strongly suggest to attend for good understanding of Modal Analysis.
