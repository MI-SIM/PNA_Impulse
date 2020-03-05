# PNA_Impulse
XPP and Matlab code to accompany the manuscript "Bifurcation analysis of an impulsive system describing Partial Nitritation and Anammox in a hybrid reactor"

## Code developed and maintained by Dr. Matthew J. Wade (matthew.wade@ncl.ac.uk, dr.matthewwade@googlemail.com)

## Files contained in this repository

Impulse.ode : An XPP (http://www.math.pitt.edu/~bard/xpp/xpp.html) file that can be used to simulate the impulsive system of ODEs and find the initial values of the state variables on the positive impulse points (period orbit). Run for t long enough to guarantee convergence to the periodic orbit. This may be confirmed by inspecting the Poincaré Map intersecting the period orbit such until sequential points are constant.

To use the following code, ensure that MatCont for maps (https://sourceforge.net/projects/matcont/files/matcontm/matcontm5p4/) is installed and added to your Matlab path.

## Main bifurcation code to generate transcritical bifurcation curves (branch curves)

branch_curves.m : Run continuation from the initial impulse points identified using XPP, find the transcritical bifurcation and follow the curve in two parameter space

pna_sys.m : The system of impulsive Ordinary Differential Equations (PN/A model), executed within branch_curves.m

## Code to generate contour maps of the performance indices. NB. The data generated by branch_curves.m can be used to plot the requisite curve over the contour map, as shown by Fig. 2 in the manuscript

performance_init.m : Initialise two parameter space to calculate performance metrics (Hydraulic Retention Time and Nitrogen removal efficiency) in the admissable operating region

ODE_performance_grid.m : Calculation of performance metrics in selected parameter space. Plots a contourmap of the resulting values

pna_sys_perf.m : Impulsive system ODE solver for performance metrics

ImpulseA.m : Calculation of the impulsive map, identical to the XPP code (Impulse.ode) but here used to solve the system until quasi steady-state to determine the system performance

RADAUsolver.m : Third-party numerical solver for stiff problems (Runge-Kutte first-order Differential Algebraic Equation solver), numerically more stable that Matlab's in-build solvers (E. Hairer and G. Wanner, University of Geneva). Available without restriction here - https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/56162/versions/1/previews/RADAU/radau.m/index.html

## Data files

*Bifurcation curve data generated using Matlab's curve fitting toolbox to identify branch points (transcritical bifurcation) for scaling performance metric plots*

fit_model1.mat : DO vs. rAMX

fit_model2.mat : DO vs fWAS

fit_model3.mat : rAMX vs fWAS





