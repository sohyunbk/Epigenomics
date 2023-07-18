#!/bin/bash
#SBATCH --job-name=FixingBarcodes        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=400gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/4_FixingBarcodes.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/4_FixingBarcodes.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-3

SampleNameList=(A619_Re3 A619_Re4 bif3_Re3 bif3_Re4)
parser.add_argument('-BAM', "--bam_file", help="Bam file to pull reads from.", required=True, dest='bam_f')
parser.add_argument('-exp_name', "--experiment_name", help="10x config file to pull scaffold names of nuclear and non-nuclear scaffolds", required=True, dest='exp')
parser.add_argument('-output_file', "--output", help="Output file to write to. If none given, output writes to stdout.", required=False, dest='o')
parser.add_argument('-threads', "--num_threads", type=int, help="Number of threads to use for processing (default: 1).", default=1)
