#!/bin/bash
#SBATCH --job-name=Try        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=60:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Try.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Try.%j.err    # Standard error log

module load SAMtools/1.16.1-GCC-11.3.0
#!/bin/bash
Path="/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/"
# Read A.txt line by line
#!/bin/bash

cd $Path

samtools view -@ 30  --output-fmt SAM -h /scratch/sb14489/3.scATAC/2.Maize_ear/3.SortedBam/3_bif3_2_Markingpcr.bam  \
 -L ZmWUS1PromoterRegions.bed \
 -o Bif3Re2_ToRef_ZmWUSPromoterRegions.sam

#samtools view -@ 30  --output-fmt SAM -h /scratch/sb14489/3.scATAC/4.Bif3Ref/3.SortedBam/3_bif3_2_Markingpcr.bam  \
# -L ZmWUS1PromoterRegions.bed \
# -o Bif3Re2_ToBif3Ref_ZmWUSPromoterRegions.sam
A00201R:581:HC7JMDSX3:3:2278:22661:11397        81      chr1    51312   3       5M1D146M     chr9    116569571       0       ATAAGAAAAAACAGGATAGCGACAATGACTAAAAATTAAGTTTGCACTTACAGGAGAAACTGATTCATCTCTTCCATGTTGTTATACAGGTAGAGCAAAGCAGTCTTCCATTCGTCATTGGTGAAACAGTATGTTGTGTTGGGTCCTACAG   FF:FFFF:FFF:,FF::F,,FF:F:FF,FF:F,FFF,FFFF,FFFFFF:FFFFFFFFFFFFF:FF:,FFFFFFFF:F:FFFFFFF:FFFFFFFFFFFFFFFFFF:FFFF:,F:FFFFFF:FFFFFFF:FF:FFFFFFFFFFFFFFFFFFFF   NM:i:5  MD:Z:5^A35A47T44C4C12   MC:Z:79M72S     AS:i:126        XS:i:123    XA:Z:chr9,-71152323,35M1I115M,5;chr9,+3512463,151M,9;    CR:Z:GGCGTTGGTTTAGAAC   CY:Z:FF:FFFF:,,,FF:,,        CB:Z:GGCGTTGGTTTAGAAG-1 BC:Z:AGAACGCC   QT:Z::F:FFFFF   GP:i:51463   MP:i:1932976260 MQ:i:0  RG:Z:1_A619_2:MissingLibrary:1:HC7JMDSX3:3
