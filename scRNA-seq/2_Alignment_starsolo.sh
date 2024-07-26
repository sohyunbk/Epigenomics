#!/bin/bash
#SBATCH --job-name=Alignment_lowTime        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=28             # Number of CPU cores per task
#SBATCH --mem=400gb                   # Job memory request
#SBATCH --time=80:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/2_Mapping.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/2_Mapping.%j.err    # Standard error log
#SBATCH --array=0-3                   # Array range

Read1=(Sohyun-wt-1_S9_R1_001.fastq.gz Sohyun-wt-2_S10_R1_001.fastq.gz Sohyun-bif3-1_S11_R1_001.fastq.gz Sohyun-bif3-2_S12_R1_001.fastq.gz)
Read2=(Sohyun-wt-1_S9_R2_001.fastq.gz Sohyun-wt-2_S10_R2_001.fastq.gz Sohyun-bif3-1_S11_R2_001.fastq.gz Sohyun-bif3-2_S12_R2_001.fastq.gz)
OutputName=(WT_Replicate1 WT_Replicate2 Bif3_Replicate1 Bif3_Replicate2)

ml Anaconda3/2020.02
source activate r_env

Path=/scratch/sb14489/4.scRNAseq/2.snRNA-seq/
Reference=/scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scRNAstarsolo/

cd "$Path"/2.Mapped_starsolo

STAR --runThreadN 28 --readFilesCommand zcat \
    --genomeDir $Reference \
    --outFileNamePrefix "${OutputName[SLURM_ARRAY_TASK_ID]}" \
    --readFilesIn $Path/1.Raw_Data/"${Read1[SLURM_ARRAY_TASK_ID]}" $Path/1.Raw_Data/"${Read2[SLURM_ARRAY_TASK_ID]}" \
    --sjdbGTFfile /scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scRNACellRanger/genes/genes.gtf \
    --soloType Droplet --soloCBwhitelist /home/sb14489/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt \
    --soloUMIlen 12 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter EmptyDrops_CR \
    --soloFeatures GeneFull \
    --soloMultiMappers PropUnique \
    --soloBarcodeReadLength 151 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 400000000000 \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
