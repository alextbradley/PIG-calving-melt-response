# Pine Island Melt Response to Calving
Code for the paper

## Drivers
MITgcm driver structures for the realistic and idealized simulations can be found in the drivers/realistic and drivers/idealized folders, respectively. The job submission scripts found within "/scripts" in each of these are appropriate for the ARCHER2 environment, where the simulations were performed.  Any input files not found there can be generated using the gendata.m scripts found in the respective folder (in "/gendata" for the idealized simulations and "/gendata_realistic" for the realistic simulations"). Note that for the idealized simulations, the topography files (which are not generated in the gendata.m script) are contained within the input folder.
