#!/bin/bash
#SBATCH --job-name=Jbrowse        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=30gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Jbrowse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Jbrowse.%j.err    # Standard error log

ml Anaconda3/5.0.1
source activate Jbrowse

## log in to Sapelo2
## Fa file and giff should have same chromosome name
NewDir=maize_ki3_WUS
cd /scratch/sb14489/0.Reference/Maize_Ki3/
/home/sb14489/Epigenomics/Jbrowse/JbrowseScripts/prepare-refseqs.pl --fasta Zm-Ki3-REFERENCE-NAM-1.0_OnlyChr_Bif3.fa --out "$NewDir"
/home/sb14489/Epigenomics/Jbrowse/JbrowseScripts/flatfile-to-json.pl --gff Zm-Ki3-REFERENCE-NAM-1.0_Zm00029ab.1_OnlyChr.gff3 --trackLabel genes --out "$NewDir"
/home/sb14489/Epigenomics/Jbrowse/JbrowseScripts/generate-names.pl -v --out "$NewDir"
#Using 1 chars for sort log names (16 sort logs) --> means Success!!

## SCP with password
sshpass -p Schmacct5$ scp  -r /scratch/sb14489/0.Reference/Maize_Ki3//maize_ki3_WUS/ schmitzlab1@heredity.genetics.uga.edu:/data01/epigenome/JBrowse/   ## It takes long~~ time like 7 hours.....
###################################################################
## log in to Jbrowse
#ssh schmitzlab1@heredity.genetics.uga.edu
#pw: Schmacct5$
#cd /data01/epigenome/JBrowse/

#vi jbrowse.conf --> then write
#[datasets."$NewDir"]
#url = ?data="$NewDir"
#name = "$NewDir"

#mkdir "$NewDir"
#cd /data01/epigenome/JBrowse/"$NewDir" --> add
#vi tracks.conf
#[general]
#dataset_id = "$NewDir"

#cd /data01/epigenome/JBrowse/
#git add jbrowse.conf
#git add "$NewDir"
#git commit -m 'added "$NewDir"'
