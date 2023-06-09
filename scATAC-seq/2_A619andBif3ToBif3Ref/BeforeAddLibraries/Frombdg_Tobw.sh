#!/bin/bash
#SBATCH --job-name=ChangeFileForJbrowse        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=10gb                   # Job memory request
#SBATCH --time=7:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/ChangeFileForJbrowse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/ChangeFileForJbrowse.%j.err    # Standard error log
