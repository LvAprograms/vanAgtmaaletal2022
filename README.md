# vanAgtmaaletal2022
Source code for models, post-processing scripts etc. 

## scripts
Collection of python scripts used to make figures for thesis/paper.
Also contains ForceAnalysis.m, which is the main Matlab visualisation script for the models. It reads .hdf5 model output files as input.

* ForceCalculator.py - class which is inherited by GPE.py. Is able to plot topography, read marker data, etc.
* GPE.py - class inheriting from ForceCalculator. Calculates gravitational potential energy difference between two columns of rocks. Not used in current manuscript
* oceantemp.py - cooling halfspace model temperature visualisation
* slabdetachment.py - plots a dataset of detachment depths versus durations.
* GradualTransition.py - Calculates necessary node x or y coordinates to implement a gradual transition from one grid spacing to one another. Could use some automation.
* plot_dragdata.py - plots output of viscous drag calculations (from matlab script ForceAnalysis.m, can't publish here unfortunately)
* plot_slabpull.py - plots output of slab pull calcuations (from matlab script ForceAnalysis.m)
* minor_plots.py - just some convenience plotting script. 

