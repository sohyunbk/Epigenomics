#!/bin/bash
## With CellRange v2
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -SampleName|--SampleName)
                SampleName=$2
                shift
                ;;
            -Ref|--Ref)
                Reference=$2
                shift
                ;;
            -InputPath|--InputPath)
                Path=$2
                shift
                ;;
            -OutfilePath|--OutfilePath)
                OutDirName=$2
                shift
                ;;
            *)
                echo "Invalid option: $1"
                exit 1
                ;;
        esac
        shift
    done

    if [[ -z $SampleName || -z $Reference || -z $Path ]]; then
        echo "Missing required arguments."
        exit 1
    fi
}

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Help message."
    exit 0
fi

parse_args "$@"

module load CellRanger-ATAC/2.0.0

cd "$OutDirName"

cellranger-atac count \
   --id="$SampleName"  \
   --reference="$Reference"  \
   --fastqs="$Path/$SampleName" --localcores=14
