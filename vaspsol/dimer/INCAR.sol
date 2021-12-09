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

 # encut = Specifies the cutoff energy for the plane wave basis set in eV.
 # nelect = Sets the number of electrons.
 # amix = specifies the linear mixing parameter.
 # ediff = specifies the global break condition for the electronic selfconsistent-loop in ev.
 # ediffg = defines the break condition for the ionic relaxation loop.
 # prec = specifies the "precision" mode.
 # istart = determines whether or not to read the wavecar file (1 = continuation job: "restart with constant energy cut-off").
 # nelm = sets the maximum number of electronic selconsistency steps.
 # isym = determines the way that vasp treats symmetry (0 = vasp does not use symmetry).
 # ispin = specifies spin polarisation (1 = non spin polarised calculations are performed).
 # sigma = specifies the width of smearing in ev.
 # ismear = determines how the partial occupancies are set for each orbital (-1 = fermi smearing).
 # algo = specifies the electronic minimisation algorithm.
 # maxmixs = pecifies the max number steps stored in broyden mixer.
 # lcharg = determines whether the charge densities (chgcar and chg) are written.
 # lwave = determines whether the wavefunctions are writen to the wavecar file at the end of a run.
 # lreal = determines whether the projection operators are evaluated in real-space or in reciprocal space.
 # nsw = sets maximum number of ionic steps.
 # ibrion = determines how the ions are updated and moved.
 # potim = sets the step width scaling for ionic relaxations.
 # iopt = choose which force based optimizer to use (2 = conjugate gradient, 3 = quick-min, 7 = fast inertial relax engine).


 # Solvent tags

 # lsol = turn solvent on or off.
 # lambda_d_k = the debeye length in angstroms.
 # eb_k = the relative permativity of the solvent.

 # Functional tags

 # gga = specifies the type of generalised-gradient-approximation one wishes to us.
 # luse_vdw = switches vdw-df functionals on and off.
 # zab_vdw = 

 # Dimer tags

 # ichain = use the dimer method.
 # ddr = the dimer seperation (twice the distance between images).
 # drotmax = maximum number of rotation steps per translation step.
 
