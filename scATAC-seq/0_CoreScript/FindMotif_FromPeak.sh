#!/bin/bash
Meme_fromPeakBedFile() {
    ##### Start
    CMD="bedtools getfasta -fi $ref -bed $bed > $bed.fa"
    eval $CMD
    CMD2="meme-chip $bed.fa -minw 4 -maxw 15 -o $Outfile"
    eval $CMD2
}

get_parser() {
    cat <<EOF
    Usage: $0 [options]

    Call Peaks for scATAC data.
    Should load modules:
    module load MEME/5.4.1-foss-2019b-Python-3.7.4
    module load BEDTools/2.30.0-GCC-10.2.0

    Options:
        -BedFile, --BedFile  BedFile
        -Ref, --Ref  Ref.FaFile
        -Outfile, --Outfile  Outfile
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -BedFile|--BedFile)
                bed=$2
                shift
                ;;
            -Ref|--Ref)
                fai=$2
                shift
                ;;
            -Outfile|--Outfile)
                Outfile=$2
                shift
                ;;
            *)
                echo "Invalid option: $1"
                exit 1
                ;;
        esac
        shift
    done

    if [[ -z $bed || -z $ref || -z $Outfile ]]; then
        echo "Missing required arguments."
        exit 1
    fi
}

parse_args "$@"
Meme_fromPeakBedFile
