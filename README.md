
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/sainag2473/ModalAnalysis_MAC_and_Modalupdating/master)



This is a study to explore Modal Analysis, Correlation and Modal Updating of FE models using Test data. In addition, Modal Analysis using Polyreference Technique for Symmetric bodies.

The project is systematically explained, can be viewed here: [Python Jupyter Notebook](https://github.com/sainag2473/MAC_and_Modalupdating/blob/master/MAC_JupyterNotebook/MACandModelCorrection.ipynb)


Edit the Jupyter notebook live in the browser. Click on the link to edit the Jupyter notebook, no python installation required, the code can be run live, any corrections of typos or wrong statements made and any addition of signal processing of the data such as using filters etc,. are welcomed for a pull request  [Find MACandModelCorrection.ipynb here to edit](https://mybinder.org/v2/gh/sainag2473/ModalAnalysis_MAC_and_Modalupdating/master)




**Product Development Process:** There are five phases in this process along the cost of change curve.

- Top two phases: **Concept and Detail Drawing**, in these phases, the CAD and FEA simulation play the major role. So, any change of product design is very inexpensive and it doesn't take much time for modifications.


- Next two phases: **Prototype and Production**, in these phases which are sort of testing phase, it becomes little bit more expensive to change the aspects of the product's design. So, it gets more expensive as moving forward in time.


- Final phase: **Field Failure**, this is the phase where product recalls occur to make a major design changes and it gets much more expensive.


Serious field failure in the final phase is because often times the biggest loss of opportunity is correlating the test data obtained in the 3rd-4th phases to the simulation model in 1st-2nd phases, this results in a FE model that would not accurately represents the real world. Because of the fact that the FE model is used in conducting different analysis such as vibration analysis, loads analysis, acoustics, durability etc., it is important prior to doing these analysis that FE model is updated using test data so that the model is developed a head of time to represent real world. Hence, we need to have the base FE model as accurate as possible as far as Mass, Stiffness, Damping, Natural Frequency, and mode shapes. 


Therefore, it is important to combine the testing and the simulation so that validating not only the physical design but also whether the assumption made for the FE simulation is correct. So, the first preference is given to validation at the test condition to use the right modeling technique in simulation to look at the variations and develop the model a head of time to see certain phenomenon that would replicate the real world.


**Repository contains details of the following**:

FE data

Testing was done in a laboratory environment using an open source experimental modal analysis software called "OpenModal" which had been validated with existing Bobcat commercial software.

### Other related projects:

[Wing Modelled using Brick Element](https://github.com/sainag2473/Project_2) 


[Wing Modelled using Beam Element](https://github.com/sainag2473/Project_1) 



## References to Modules and Books Used:

[WFEM](https://github.com/josephcslater/WFEM) : by Dr. Joseph Slater

[OpenModal](https://github.com/openmodal/OpenModal) : created as an Opensource software written in 
Python for Vibration Testing with excellent capabilities. 

[vibrationtesting](https://github.com/Vibration-Testing/vibrationtesting) : by Dr. Joseph Slater

[Vibration Toolbox]( https://github.com/vibrationtoolbox/vibration_toolbox.git) : by Dr. Joseph Slater

Slater, J. C., Vibration Testing, with Modal Analysis and Health Monitoring, John Wiley & Sons, Ltd, 2002
  
Finite Element Model Updating in Structural Dynamics by M.I. Friswell and J.E. Mottershead

Concepts and Applications of Finite Element Analysis, Fourth edition by Robert D. Cook, David S. Malkus, Michael E. Plesha and Robert J. Witt

[Siemens Corp. Seminar](https://community.plm.automation.siemens.com/t5/Event-Collateral/Modal-Analysis-and-FEA-Test-Correlation/ta-p/435288): Presenters were Mr Scott Beebe and Ms Jacklyn Kinsler

[Bruel & Kjaer Sound and Vibration Conference](https://bksv.com/en/Training/training-courses/UnitedStates/UnitedStates/Live-Courses/Modal-Analysis/05-15-2018): Some Pictures and Bulletin points were taken directly from Modal Analysis Conference conducted by Bruel and Kjaer, Chicago, May2018.
Mr Robert Trepanier is the presenters who has more than 30 years of experience, I strongly suggest to attend for good understanding of Modal Analysis.


## Copyright

 This book is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). you can freely use this notebook for educational/academic purposes.
