#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --mem=256G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64


module load StdEnv/2023
module load smrtlink/25.3
module load samtools



./rnaquant.sh 64
