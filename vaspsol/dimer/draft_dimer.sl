#!/bin/bash -e

#SBATCH --job-name        graph_Hey_diff_dimer_hugepages_potstep_1                                
#SBATCH --account         uoo02693
#SBATCH --time            24:00:00
#SBATCH --nodes           1
#SBATCH --ntasks-per-node 40               # One CPU per task is assumed
#SBATCH --mem-per-cpu     1000
#SBATCH --hint            nomultithread
#SBATCH --output          single_pot_step.%j.out # Include the job ID in the names of
#SBATCH --error           single_pot_step.%j.err # the output and error files
#SBATCH --mail-user=warci626@gmail.com
#SBATCH --mail-type=END,FAIL,REQUEUE

ulimit -s unlimited
module load gcc/6.1.0
module load VASP/5.4.4-CrayIntel-19.04b-VTST-sol-hugepages
module switch craype-hugepages2M craype-hugepages8M

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
nel=`grep NELECT INCAR | awk '{print $3}'`
echo 'NELECT    FermiShift    FermiE    SystPot' 	
echo   $nel  $chift  $fermi  $udiff

if (( $(echo "sqrt($udiff^2) < $thr" | bc -l) )); then
    echo 'Target potential reached'
fi

nelect=`echo "$nel + $factor * $udiff" | bc -l`

echo $nelect

#  sed -i "/NELECT/c\ NELECT = $nelect" INCAR.vac
#  sed -i "/NELECT/c\ NELECT = $nelect" INCAR.sol

rm WAVECAR
