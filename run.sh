#!/bin/bash

#SBATCH --job-name=NvP53
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

python ./run_NvP53.py ../dataset/Total-20220906/Dataset.pmat.pb.20220906.h5