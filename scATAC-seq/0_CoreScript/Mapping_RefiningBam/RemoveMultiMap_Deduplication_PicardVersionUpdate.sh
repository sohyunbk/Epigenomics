#!/bin/bash

### This step is after mapping
## Should load this packages
## It's for the output of cellranger v2

module load Anaconda3/2022.10
source activate /home/sb14489/.conda/envs/r_env
module load picard/2.27.5-Java-15
module load  SAMtools/1.10-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0
# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --path)
            Path="$2"
            shift
            shift
            ;;
        --MappedDir)
            MappedDir="$2"
            shift
            shift
            ;;
        --OGSampleName)
            OGSampleName="$2"
            shift
            shift
            ;;
        --NewSampleName_forBam)
            NewSampleName_forBam="$2"
            shift
            shift
            ;;
        *)
            echo "Unknown option: $key"
            exit 1
            ;;
    esac
done

# Function 1: Remove Multimap
remove_multimap() {
    mkdir -p 3.SortedBam

    samtools view -@ 24 -h -f 3 -q 10 "$Path"/"$MappedDir"/"$OGSampleName"/outs/possorted_bam.bam |
    grep -v -e 'XA:Z:' -e 'SA:Z:' |
    samtools view -@ 24 -bS - > "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted.bam
    ### -v -e : exclude the lines
    # XA:Z: SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ Other canonical alignments in a chimeric alignment,
    # XA: Alternative hits; format: (chr,pos,CIGAR,NM;)*
    # SA: Z, not sure what it is but, it almost always coincides with the 256 flag = not primary alignment

    ##########################################################################################
    ## Check and change Header for PICARD and Cellrangerv2
    ##########################################################################################
    samtools view -H "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted.bam > "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted_Header.sam
    sed -i '1d' "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted_Header.sam
    sed -i '1i @HD\tVN:1.5\tSO:coordinate' "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted_Header.sam

    samtools reheader "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted_Header.sam  "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted.bam \
    > "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted_HF.bam

    samtools view -h "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted_HF.bam | \
    awk 'BEGIN {OFS="\t"} { if ($0 ~ /^@/) {print $0} else { gsub(/CB:Z:[^-\t]+-1/, "CB:Z:\\1", $0); print $0 } }' | \
    samtools view -b -o "$Path"/3.SortedBam/"$NewSampleName_forBam"_Sorted_HF_FixingCBBarcode.bam

}

# Function 2: Deduplication
deduplication() {
    cd "$Path"/3.SortedBam
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 MAX_RECORDS_IN_RAM=1500000\
        REMOVE_DUPLICATES=true METRICS_FILE="./$NewSampleName_forBam"_dups_Markingpcr.txt \
        I="./$NewSampleName_forBam"_Sorted_HF_FixingCBBarcode.bam \
        O="./$NewSampleName_forBam"_Rmpcr.bam \
        BARCODE_TAG=CB \
        ASSUME_SORT_ORDER=coordinate \
        USE_JDK_DEFLATER=true
}

# Call the functions
remove_multimap
deduplication
