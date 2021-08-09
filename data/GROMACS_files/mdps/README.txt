This directory contains the .mdp files used in GROMACS simulations.
There are nine mdp files:

  ** for relaxing molecules in the initial (gas phase) configuration and 
     compressing into a (poorly equilibrated) melt **

  1. g_em.mdp        : energy minimization
  2. g_nvt_short.mdp : 0.1 fs timestep NVT
  3. g_nvt.mdp       : 1.0 fs timestep NVT
  4. g_compress.mdp  : 100 bar NPT

  
  ** for compressing the initial melt to a high density, eliminating any voids
     from incomplete compression in step 4 **

  5. c1_nvt.mdp : high temperature (1000 K) NVT
  6. c2_nvt.mdp : room temperature (300 K) NVT
  7. c3_npt.mdp : room temperature, high pressure (1000 bar)  NPT
  
  ** burn-in period to de-equilibrate the melt from its post-compression state.
     switches barostat to Parinello-Rahman for valid sampling **

  8. burn_in.mdp


  ** production annealing simulations **

  9. anneal_schedule.mdp

