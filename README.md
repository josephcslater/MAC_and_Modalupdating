# Model Assurance Criterion and Model correction or updating: 

This is an Independent Research study to explore Modal Updating of FE models using Test data. 

Purpose of bridging the gap between NVH and CAE: Product development process have been modifying every now and then, at present most of the companies don't use simulation unless they first establish trust. So, first preference is given to validation at the test condition to use the right modeling techniques in simulation so that the variations can be looked at and the model could be used for development a head of time to see certain phenomenon that would replicate real world. 

Repository contains details of the following:
FE models using 6 d.o.f/node Beam element and 8 noded Brick element were coded in MATLAB which had been validated with Ansys Apdl. 

Testing was done in a laboratory environment using an open source experimental modal analysis software called "OpenModal" which had been validated with existing Bobcat commercial software.

The project is systematically explained in the [Python jupyter notebook](https://github.com/sainag2473/MAC_and_Modalupdating/blob/master/MAC_JupyterNotebook/MACandModelCorrection.ipynb)

## References to Modules and Books Used:

[WFEM](https://github.com/josephcslater/WFEM) : by Dr. Joseph Slater

[OpenModal](https://github.com/openmodal/OpenModal) : created as an Opensource software written in 
Python for Vibration Testing with excellent capabilities. 

[vibrationtesting](https://github.com/Vibration-Testing/vibrationtesting) : by Dr. Joseph Slater

[Vibration Toolbox]( https://github.com/vibrationtoolbox/vibration_toolbox.git) : by Dr. Joseph Slater
  
Finite Element Model Updating in Structural Dynamics by M.I. Friswell and J.E. Mottershead

Concepts and Applications of Finite Element Analysis, Fourth edition by Robert D. Cook, David S. Malkus, Michael E. Plesha and Robert J. Witt

[Siemens Corp. MAC](https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Modal-Assurance-Criterion-MAC/ta-p/368008)
