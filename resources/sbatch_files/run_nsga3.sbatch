#!/bin/bash 
# 
#SBATCH --job-name=run_simulation
# 
#SBATCH --time=72:00:00
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1

#SBATCH --partition serc

#SBATCH --mem=3G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail 
#SBATCH --mail-user=jshs@stanford.edu

# execute script
cd /scratch/users/$USER/swat-pytools/tests/
eval "$(conda shell.bash hook)"
conda activate swatpy
./run_optimization.py
