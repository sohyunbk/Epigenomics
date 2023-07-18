##Transfer python script to #!/usr/bin/env bash
#!/bin/bash

function from_bdgfile_to_bwfile() {
    local BdgFile=$1
    local OutFileName=$2
    local Fai=$3

    Cmd_sort="bedSort ${BdgFile} ${BdgFile}_Sorted"
    Cmd="bedGraphToBigWig ${BdgFile}_Sorted ${Fai} ${OutFileName}.bw"

    $Cmd_sort
    $Cmd
}

function from_bedfile_to_dirforTrack() {
    local BedFile=$1
    local OutFileName=$2
    local Path=$(dirname $OutFileName)
    local Prefix=$(basename $OutFileName)

    Cmd="/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed ${BedFile} --trackLabel ${Prefix} --out ${Path}"
    $Cmd
}

function make_bed_from_samfile() {
    local Samfile=$1
    local ReadLength=$2
    local OutFileName=$3

    Outfile=$OutFileName
    Infile=$Samfile

    while IFS=$'\t' read -r -a List; do
        nFragment=${List[8]}
        nFragment=${nFragment/#-/}
        if (( nFragment < ReadLength )); then
            nlength=$nFragment
        else
            nlength=$ReadLength
        fi
        nStart=${List[3]}
        sChr=${List[2]}
        echo -e "$sChr\t$nStart\t$((nStart + nlength))" >> $Outfile
    done < $Infile
}

function from_bam_to_bwfile() {
    local Bamfile=$1
    local Fai=$2

    samtools index -@ 24 $Bamfile
    python /home/bth29393/jbscripts/file_to_bigwig_pe.py $Fai $Bamfile
    BedName="${Bamfile%.bam}.bed"
    bedtools bamtobed -i $Bamfile > $BedName
    bedtools genomecov -i "$BedName" -split -bg -g "$Fai"  > "${Bamfile%.bam}.bg"
    wigToBigWig "${Bamfile%.bam}.bg" "$Fai"  "${Bamfile%.bam}.bw"
}

Step="${Step:-}"
OutputName="${OutputName:-}"
BdgFile="${bdgFile:-}"
Fai="${Fai:-}"
BedFile="${bed:-}"
Samfile="${sam:-}"
Bamfile="${bam:-}"
ReadLength="${readlength:-}"

if [[ $Step == "bdgTobw" ]]; then
    from_bdgfile_to_bwfile $BdgFile $OutputName $Fai
elif [[ $Step == "BedToTrack" ]]; then
    from_bedfile_to_dirforTrack $BedFile $OutputName
elif [[ $Step == "SamToBed" ]]; then
    make_bed_from_samfile $Samfile $ReadLength $OutputName
elif [[ $Step == "BamTobw" ]]; then
    from_bam_to_bwfile $Bamfile $Fai
fi

#bash script.sh -Step <Step> -OutputName <OutputName> -bdgFile <BdgFile> -Fai <Fai> -bed <BedFile> -sam <Samfile> -bam <Bamfile> -readlength <ReadLength>
