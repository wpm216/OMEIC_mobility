# OMEIC_mobility
Hole mobility in Organic Mixed Electronic and Ionic Conductors (OMEICS),
accompanying 
'Efficient electronic transport through the insulator-rich phase in PEDOT:PSS blends'
by S.T. Keene, W. Michaels, et al.

## Contents:
#### `/data/`
- `/data/distance_matrices`: 
    inter-cluster distances of PEDOT, calculated with atomistic MD.
    See SI Section S1 for generation details.
    and SI Sections S2-S3 for usage.
- `/data/island_MC_sims_varying_dilution`: 
    configurations of PEDOT-rich grains generated by MC simulation.
    See SI Section S4 for generation details
    and SI Section S3 for usage.
- `/data/allrates.mat`: 
    Matlab datafile containing all data used in the publication.
    Array and variable names are used in `/src/*` scripts.

#### `/src/`
Contains analysis code to generate Figure 4e, Figure 4f,
Extended Data Fig. 6, and Extended Data Fig. 7.
- Files beginning with `production`, e.g. `/src/productionAtomistic.m`
generate steady-state charge-transfer rates between pairs of PEDOT-rich grains.
(see SI Section S3).
- Files beginning with `plot`, e.g. `/src/plotFigure4e`,
calculate film mobility vs. PEDOT to PSS weight ratio
(see SI Section S3).
- Other files contain accessory functions.

