#!/bin/bash
module load BEDTools/2.30.0-GCC-10.2.0
# vars
name=$1
peaks=$( ls *_summits.passing_FDR_Edited )

# iterate
for i in ${peaks[@]}; do

	echo "$i"

	id=$( echo "$i" | cut -d'.' -f1-2 )

	perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-2_normalize_score.pl $i \
	| perl -ne 'chomp;my@col=split("\t",$_);
		$col[1] = $col[1] - 250;
		$col[2] = $col[2] + 250;
		if($col[1] >= 0){
			print "$col[0]\t$col[1]\t$col[2]\t$col[3]\t$col[0]_$col[1]_$col[2]_$col[3]\n";}' - > $id.temp
done

# merge peaks
cat *.temp | sort -k1,1 -k2,2n - > $1.temp2

# actually merge
bedtools merge -i $1.temp2 -c 5 -o collapse | perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-2_selectNonOverlapping.pl - > $1.unique500bpPeaks.bed
rm *.temp *.temp2
