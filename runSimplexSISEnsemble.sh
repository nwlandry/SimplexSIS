#!/bin/bash

#SBATCH --qos=blanca-appm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=nila0959@colorado.edu
#SBATCH --job-name=SimplexSIS
#SBATCH --error=slurm.%j.err
#SBATCH --output=slurm.%j.out

module purge
module load slurm/blanca
source /curc/sw/anaconda3/2019.03/bin/activate
conda activate idp

python runModelEnsemble.py
