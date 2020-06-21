#!/bin/bash

#SBATCH --qos=blanca-appm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=143:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=nila0959@colorado.edu
#SBATCH --job-name=SimplexSIS
#SBATCH --error=slurm.%j.err
#SBATCH --output=slurm.%j.out

module purge
module load slurm/blanca
source /curc/sw/anaconda3/2019.03/bin/activate
conda activate idp

python runModelInParallel.py 10000 "uniform" 50 150 True 100 4
