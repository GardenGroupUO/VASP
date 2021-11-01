 # Dimer INCAR.sol
 
 # Standard VASP Tags

 ENCUT = 500.000000
 NELECT = 1
 AMIX = 0.050000
 EDIFF = 1.00e-07
 EDIFFG = -0.05
 PREC = Normal
 ISTART = 1
 NELM = 120
 ISYM = 0
 ISPIN = 1
 SIGMA = 0.10000
 ISMEAR = -1
 ALGO = Fast
 MAXMIX = 80
 LCHARG = .FALSE.
 LWAVE = .FALSE.
 LREAL = Auto
 NSW = 10000
 IBRION = 3
 POTIM = 0
 IOPT = 2 

# Solvent Tags

 LSOL = .TRUE.
 LAMBDA_D_K = 3.000000
 EB_K = 78.400000

 
# Functional Tags

 GGA = BF
 LUSE_VDW = .TRUE.
 Zab_VDW = -1.8867 
 
 # Dimer Tags

 ICHAIN = 2
 DdR = 1e-2 
 DRotMax = 6



 # Standard VASP tags

 # ENCUT = Specifies the cutoff energy for the plane wave basis set in eV.
 # NELECT = Sets the number of electrons.
 # AMIX = Specifies the linear mixing parameter.
 # EDIFF = Specifies the global break condition for the electronic Selfconsistent-loop in eV.
 # EDIFFG = Defines the break condition for the ionic relaxation loop.
 # PREC = Specifies the "precision" mode.
 # ISTART = Determines whether or not to read the WAVECAR file (1 = Continuation job: "restart with constant energy cut-off").
 # NELM = Sets the maximum number of electronic selconsistency steps.
 # ISYM = Determines the way that VASP treats Symmetry (0 = VASP does not use symmetry).
 # ISPIN = Specifies spin polarisation (1 = non spin polarised calculations are performed).
 # SIGMA = Specifies the width of smearing in eV.
 # ISMEAR = Determines how the partial occupancies are set for each orbital (-1 = Fermi smearing).
 # ALGO = Specifies the electronic minimisation algorithm.
 # MAXMIXS = pecifies the max number steps stored in Broyden mixer.
 # LCHARG = Determines whether the charge densities (CHGCAR and CHG) are written.
 # LWAVE = Determines whether the Wavefunctions are writen to the WAVECAR file at the end of a run.
 # LREAL = Determines whether the projection operators are evaluated in real-space or in reciprocal space.
 # NSW = Sets maximum number of ionic steps.
 # IBRION = Determines how the ions are updated and moved.
 # POTIM = Sets the step width scaling for ionic relaxations.
 # IOPT = Choose which force based optimizer to use (2 = Conjugate Gradient, 3 = Quick-Min, 7 = Fast Inertial Relax Engine).


 # Solvent tags

 # LSOL = Turn solvent on or off.
 # LAMBDA_D_K = The Debeye length in Angstroms.
 # EB_K = The relative permativity of the solvent.

 # Functional tags

 # GGA = Specifies the type of generalised-gradient-approximation one wishes to us.
 # LUSE_VDW = Switches vDW-DF functionals on and off.
 # Zab_VDW = 

 # Dimer tags

 # ICHAIN = Use the dimer method.
 # DdR = The dimer seperation (twice the distance between images).
 # DRotMax = Maximum number of rotation steps per translation step.
 