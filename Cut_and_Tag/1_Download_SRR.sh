#!/bin/bash
#SBATCH --job-name=SRAToolKit        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=30gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/SRAToolKit.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/SRAToolKit.%j.err    # Standard error log

module load SRA-Toolkit/3.0.3-gompi-2022a

while IFS= read -r line; do
    srr_id=$(echo $line | awk '{print $1}')

    # Skip the header or any empty lines
    if [[ $srr_id == SRR* ]]; then
        echo "Downloading $srr_id..."
        fasterq-dump $srr_id --split-files --outdir /scratch/sb14489/7.DAPorChIP/CUTandTAG/1.RawData/
    fi
done < /scratch/sb14489/7.DAPorChIP/CUTandTAG/1.RawData/Cut_And_Tag_SRR.txt

echo "Downloads complete."
