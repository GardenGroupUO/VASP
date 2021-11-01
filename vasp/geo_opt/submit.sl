#!/bin/bash -e
#SBATCH -J edge_Mo
#SBATCH -A uoo03025         # Project Account
#SBATCH --time=23:59:00     # Walltime
#SBATCH --ntasks=12 #number of nodes * tasks(12)
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-user=warci626@gmail.com
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --hint=nomultithread  #no hyperthreading


module load gcc/8.3.0
module load VASP/5.4.4-CrayIntel-18.08-VTST-sol
srun vasp_std
