#!/bin/bash -e

#SBATCH --job-name        3w_Au_spaced_S_Hads_VASPsol                                 
#SBATCH --account         uoo02693
#SBATCH --time            24:00:00
#SBATCH --nodes           1
#SBATCH --ntasks-per-node 40               # One CPU per task is assumed
#SBATCH --mem-per-cpu     1000
#SBATCH --hint            nomultithread
#SBATCH --output          MyVASPJob.%j.out # Include the job ID in the names of
#SBATCH --error           MyVASPJob.%j.err # the output and error files
#SBATCH --mail-user=warci626@gmail.com
#SBATCH --mail-type=END,FAIL,REQUEUE

ulimit -s unlimited
module load gcc/6.1.0
module load VASP/5.4.4-CrayIntel-19.04b-VTST-sol
#SUBMITDIR=/nesi/nobackup/uoo02693/$SLURM_JOB_ID
#mkdir $SUBMITDIR
#cp -pr ./* $SUBMITDIR/
#CURRENT_DIR=$PWD
#cd $SUBMITDIR

x=0
for i in {0..99}
do
        if [ ! -d results-$i ]
        then
                mkdir results-"$i"
                x=$i
                break
        else
                continue
        fi
done

# This script will run several constant-charge
# dimer runs in a row, with intermittent updates
# on the number of electrons so that the potential
# at the saddle point converges to the chosen
# target potential.

# stopping threshold in V based on |U-$utarget|
thr=0.002  
# target potential in V vs $phiref
utarget=0.0 

#ncores=32  # number of CPUs
# absolute reference potential in V vs vacuum
phiref=4.43  
# determines how to calculate the next number
factor=1.0  
# of electrons NELECT based on the current number
# and the deviation from the target potential:
# NELECT(i+1) = NELECT(i) + factor * |U-$utarget|
# This factor is hence related to the capacitance and
# the surface area. Simply taking factor=1 worked well enough
# for Cu and Pt slabs with 12 metal atoms at the surface.

# In addition to the KPOINTS, POSCAR, POTCAR and MODECAR
# files, the script requires 2 INCAR files:
# INCAR.vac: for a single-point run without solvent to
#            write the wave functions
# INCAR.sol: for a dimer run with the solvent
# Note that, for this script to work, both INCAR files
# need to explicitly provide the (same) initial guess
# for the number of electrons (the NELECT tag).
for i in {0..99}; do
  cp INCAR.vac INCAR

#  mpirun -np $ncores vasp_std 1>out 2>err

  srun vasp_std 1>out 2>err

  mv OSZICAR OSZICAR.vac
  mv OUTCAR OUTCAR.vac
  cp INCAR.sol INCAR

  srun vasp_std 1>out 2>err
#  mpirun -np $ncores vasp_std 1>>out 2>>err

  chift=`grep -i FERMI_SHIFT out | tail -n 1 | cut -d= -f 2`
  fermi=`grep -i fermi OUTCAR | tail -1 | awk '{print $3}'`
  udiff=`echo "-1*($phiref + $chift + $fermi) - $utarget" | bc -l`
  pot=`echo "$phiref + $chift + $fermi"`
  nel=`grep NELECT INCAR | awk '{print $3}'`
  echo 'Iteration    NELECT    FermiShift    FermiE    PotentialDiff'
  echo $i  $nel  $chift  $fermi  $udiff

  if (( $(echo "sqrt($udiff^2) < $thr" | bc -l) )); then
    break
  fi

  nelect=`echo "$nel + $factor * $udiff" | bc -l`
  sed -i "/NELECT/c\ NELECT = $nelect" INCAR.vac
  sed -i "/NELECT/c\ NELECT = $nelect" INCAR.sol
  
  mkdir results-$x/potstep-$i
  mv INCAR results-$x/potstep-$i/INCAR-$i
  mv POSCAR results-$x/potstep-$i/POSCAR-$i
  mv OUTCAR results-$x/potstep-$i/OUTCAR-$i
  mv out results-$x/potstep-$i/out-$i
  mv err results-$x/potstep-$i/err-$i
  mv OSZICAR results-$x/potstep-$i/
  mv OSZICAR.vac results-$x/potstep-$i/
  mv PCDAT results-$x/potstep-$i/
  mv REPORT results-$x/potstep-$i/
  mv XDATCAR results-$x/potstep-$i/
  mv vasprun.xml results-$x/potstep-$i/
  mv EIGENVAL results-$x/potstep-$i/
  mv IBZKPT results-$x/potstep-$i/
  mv CHG results-$x/potstep-$i/
  mv CHGCAR results-$x/potstep-$i/
  mv DOSCAR results-$x/potstep-$i/
  mv OUTCAR.vac results-$x/potstep-$i/
  mv CONTCAR CONTCAR-$i


  cp CONTCAR-$i results-$x/potstep-$i/CONTCAR-$i


  mv CONTCAR-$i POSCAR
  rm WAVECAR
done
