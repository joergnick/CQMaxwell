
The following scripts lead to the reconstruction of the images in the preprint "Time-dependent electromagnetic scattering from thin layers". 

RKRefErrorDatadelta01.py
RKRefErrorDatadelta10.py

These two scripts (runtime each around 3-4 days) generate the data for the convergence plots. Then, the Matlab Script 

/PlottingScripts/ErrorPlotsMain.m

generates the convergence plots.

The data related to the other figures is generated from 

FramesAndConditions.py . (Runtime around 3-4 days as well.)


Then, the remaining figures are generated by

/PlottingScripts/FramesPlot.m
/PlottingScripts/ConditionContourPlot.m.

###################################################################################################

Possible memory issue with FramesAndConditions.py :

Got a Segmentation fault at first, used 

ulimit -s 20000,

to increase the stack allocated to the python process. Originally, the output of ulimit -s was 8192. Might be advisable as a precaution.

Used Versions of some relevant packages : 

Bempp version : 3.3.5
Scipy version : 0.17.0
Numpy version : 1.11.0
Matplotlib version : 1.5.1

