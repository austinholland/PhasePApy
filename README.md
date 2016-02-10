
========================================================================================
			PhasePApy Version 1.1
========================================================================================

			
			### ABOUT PhasePApy ###

The PhasePApy is a Pure Python program package developed by Chen Chen (c.chen@ou.edu; 
c.chen8684@gmail.com) under the guidance of Austin A. Holland (austin.holland@ou.edu; 
aaholland@usgs.gov).

Version 1.0: 2014	--	FBpicker; 1D Associator
Version 1.1: 2016	--	FBpicker, AICDpicker, KTpicker; 1D and 3D Associator

The PhasePApy is a Seismic Phase Picker and Associator program package, which is 
designed to identify and associate picks to the seismic events or microseismic events. 
This work has currently been submitted for peer review. Please reference the paper in 
review if you are using this code. Once the paper is published we will include the 
full reference.

The first version 1.0 has been used by Oklahoma Geological Survey (OGS) for near real
time earthquake monitoring since 2014. Here, we only provide the version 1.1 to public,
because we think this one is more user-friendly and easy to manipulate. You might need 
to modify the code a bit for your requirements. Also, we still think the code has some 
space for further clean-up and improvement, so please feel free to modify and update it. 

If you have certain ideas to improve this package, please let us know.
If you are going to use this package for any commercial purpose, please let us know.


			### COMPONENTS OF PhasePApy ###

The PhasePApy package composes of two sub-packages: PhasePicker and Associator. These 
two sub-packages can work jointly or separately, depending on users’ requirements. 
For version 1.1, the PhasePicker contains three pickers: FBpicker, AICDpicker, and KTpicker;
the Associator includes the 1D and 3D Associator modules. The PhasePicker and Associator is 
designed to identify phase onsets and associate them to best fitting events. P- and S-phases
can be associate in this version, but the algorithm can be accommodated for any distance.


			### EXAMPLES ###

In the test directory, there are two examples to show how to tune the parameters and present 
the performance with all the processing units. The first example is the 6-minute data and the
second one is a one-day dataset. The README.md file is included in each example. 


			### REFERENCES ###

Chen,C., and A. A. Holland, (201X). PhasePApy: A Robust Pure Python Package for Automatic 
Identification of Seismic Phases, XXX XXX XXX, XX(XX), XXX-XXX.
			

			### REQUIRED PYTHON LIBRARIES ###	

The PhasePApy package relies on open libraries: Obspy, Numpy, and Scipy. if you are going to 
plot results with plotting scripts in this package, you need matplotlib as well. Otherwise, 
the users can visualize it based on their own methods.	




