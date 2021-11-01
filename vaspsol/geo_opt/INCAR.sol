 # VASPsol geo-opt INCAR.sol

 # Standard VASP Tags

 ENCUT = 500.000000
 NELECT = 1
 AMIX = 0.050000
 EDIFF = 1.00e-05
 EDIFFG = -0.03
 PREC = Normal
 ISTART = 1
 NELM = 120
 ISYM = 0
 ISPIN = 1
 SIGMA = 0.100000
 ISMEAR = -1
 ALGO = Fast
 MAXMIX = 80
 LCHARG = .FALSE.
 LWAVE = .FALSE.
 LREAL = Auto
 NSW = 1000
 IBRION = 2
 POTIM = 0.2
 LPLANE = .TRUE.
 ISIF = 0

 # Solvent Tags

 LSOL = .TRUE.
 LAMBDA_D_K = 3.000000
 EB_K = 78.400000

 # Functional Tags

 GGA = BF
 LUSE_VDW=.TRUE.
 Zab_VDW=-1.8867


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
 # LPLANE = switches on the plane-wise data distribution in real space.
 # ISIF = determines whether the stress tensor is calculated and which principal degrees-of-freedom are allowed to change in relaxation (0 = off, only forces are calculated). 
 
 # Solvent tags

 # LSOL = Turn solvent on or off.
 # LAMBDA_D_K = The Debeye length in Angstroms.
 # EB_K = The relative permativity of the solvent.

 # Functional tags

 # GGA = Specifies the type of generalised-gradient-approximation one wishes to us.
 # LUSE_VDW = Switches vDW-DF functionals on and off.
 # Zab_VDW = 

