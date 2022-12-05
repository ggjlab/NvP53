#!/bin/bash

#SBATCH --job-name=NvP53
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --output=%J.out
#SBATCH --error=%J.err

python ./prepare_pretrain.py