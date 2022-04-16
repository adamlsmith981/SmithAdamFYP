# FYP
 Source code used to generate results and plots for the FYP

### The if lucy peggs is reading this section 
Hi Lucy, Adam from the last commit to this repo here. Here is a list of what to look at and what to avoid because not everything is complete/it didnt work. 
1. MoS2, MoS2-bulk (these have lovely demonstrations on the energy bnd diagram of what a semi looks like)
2. Ta2C and Ta2CH2 (you probably have cos they are trivial)
3. Ta2CO and Ta2CO2 yer probably dont have. 
4. Keep an eye on a latch ditch effort to get both Ta2CF2 and Ta4C3 working but I would be reluctantly hopefully.
5. Just no one cares about Ta3C2H or Ta3C2OH also they dont work

## Structure of the project
In this project you will find code that was used for each material, consiting of: 
1. VC-relax
    * Including jmol display of the crystal
2. Convergence Testing
    * ecutwfc (Kinetic energy cutoff for wavefunctions, in Ry units)
    * K-points 
    * Lattice Constant
3. Kpath determined using Xcrysden
4. Bands Calculation
5. Density of States (DOS) Calculation

## Pseudopotenitals
The pseudopotentials used were the PBESOL (PBE functrional revised for solids) and takes into consideration the second order expansion for the exchange energy, which a PBE potential does not. These are Pseudopotential type: PAW. Non linear core correction. Scalar relativistic. 

Please refer to my final report where you will find more details on the computational methods. 