# The influence of Pine Island Ice Shelf calving on basal melting
This repository contains code to drive the MITgcm simulations and produce the figures contained in "The influence of Pine Island Ice Shelf calving on basal melting" by Bradley et al. (2022), DOI: 10.1029/2022JC018621. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022JC018621

## Figures
Code to produce the figures in contained in the /make_figures folder. Each figure corresponds to an individual file (for example, the file to produce figure 10 is named make_figure10.m). Figures are produced in MATLAB. NB: many of the script rely on simulation output which is not hosted in this repository because of space restrictions; to obtain a copy of these results, please contact the Alex Bradley (aleey@bas.ac.uk).

## Drivers
MITgcm driver structures for the realistic and idealized simulations can be found in the drivers/realistic and drivers/idealized folders, respectively. The job submission scripts found within "/scripts" in each of these are appropriate for the ARCHER2 environment, where the simulations were performed.  Any input files not found there can be generated using the gendata.m scripts found in the respective folder (in "/gendata" for the idealized simulations and "/gendata_realistic" for the realistic simulations"). Note that for the idealized simulations, the topography files (which are not generated in the gendata.m script) are contained within the input folder.
