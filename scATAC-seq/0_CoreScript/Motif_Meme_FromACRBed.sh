#!/bin/bash
## Load module!!
# module load MEME/5.5.0-gompi-2021b
# module load BEDTools/2.30.0-GCC-11.3.0

#Infile_Bed="/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC.FDR0.05.Bed"
#Fafile="/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa"
#MemeMotifDB="/scratch/sb14489/3.scATAC/0.Data/Plant_Motif_PWM/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt"
#OutfilePathName="/scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/IM_OC_dACR_JASPARMotif_Xstream"
#AllPeakMinusForControl="/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/Cat_A619_Bif3_500bpPeak.bed"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --infile_Bed) Infile_Bed="$2"; shift ;;
        --Fa) Fafile="$2"; shift ;;
        --MemeMotifDB) MemeMotifDB="$2"; shift ;;
        --ControlFA) ControlFA="$2"; shift ;;
        --OutfilePath) OutfilePathName="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

## 1) prepare fa
Infile_FA="$Infile_Bed".fa

if [ -e "$Infile_FA" ]; then
    echo "$Infile_FA exists."
else
    bedtools getfasta -fo "$Infile_FA" -fi "$Fafile" -bed "$Infile_Bed"
fi

#### 2)XSTREAM
xstreme  --p "$Infile_FA" --m "$MemeMotifDB" --o "$OutfilePathName" --n "$ControlFA"

## *** denovo motif search
#if [ -z "$MemeMotifDB" ]; then
#    xstreme
## ***** known motif search
#else
#    fimo --o "$OutfilePathName" "$MemeMotifDB" "$Infile_FA"
#fi
