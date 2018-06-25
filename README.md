# Inhibitor_Model
Matlab files associated with the paper: 'Predictive Modelling of a Novel Anti-adhesion Therapy to Combat Bacterial Colonisation of Burn Wounds'

https://zenodo.org/badge/137582957.svg

This folder contains the following files:

1. Command_Line_Inhibitor_Simulations_and_Plots.m: Open this file to run simulations and plot the figures. On line 36 of this code you can choose the parameter set (between 1 and 12) for which the equations are solved.

The remaining files are called by the command line when running simulations, these are:
 
 2. Fitted_Parameters: a data file containing the 13 fitted parameter values for parameter sets 1-12.
 3. No_Treatment_ODE_Fn.m: the ODE function for the untreated case with clearance in the first 24 hours.
 4. No_Treatment_Const_Deb_ODE_Fn.m: the ODE function for the untreated case with constant debridement.
 5. Inhib_Dose_ODE_Fn.m: the ODE function for when inhibitors are used, including clearance in the first 24 hours.
 6. Inhib_Dose_ODE_Fn_NC.m: the ODE function for when inhibitors are used without clearance - this is used for repeat doses in the regular inhibitor dosing scenario.
 7. Inhib_Dose_Const_Deb_ODE_Fn.m: the ODE function for when inhibitors are used with constant debridement.
