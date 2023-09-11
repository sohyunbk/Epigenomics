#!/bin/bash
#Input should look like :
#chr1	17067	17567	early_vascular_parenchyma
#chr1	17067	17567	guard_mother_cell
#chr1	17067	17567	hypodermal_sclerenchyma
#chr1	17067	17567	leaf_primordia
#chr1	17067	17567	mesophyll_precursors
#chr1	17067	17567	phloem_parenchyma
#chr1	17067	17567	procambial_meristem
#chr1	17067	17567	subsidiary_cells

## More info here: https://github.com/FunctionLab/selene/blob/master/tutorials/getting_started_with_selene/getting_started_with_selene.ipynb

## To use part of the sample:
cd /scratch/sb14489/8.ML_ACR/1.InputBed
cat phloem_SE_procambial_precursors.accessible_ACRs.bed.500 mesophyll_precursors.accessible_ACRs.bed.500 \
 axillary_meristem.accessible_ACRs.bed.500 guard_cell.accessible_ACRs.bed.500 mesophyll.accessible_ACRs.bed.500 > \
 Seedling_FiveClasses.500.bed
input_file=Seedling_FiveClasses.500.bed

conda activate /home/sb14489/miniconda3/envs/pytorch

sort -k1V -k2n -k3n "$input_file" > "$input_file".sorted

bgzip -c "$input_file".sorted > "$input_file".sorted.gz

tabix -p bed "$input_file".sorted.gz

cut -f 4 "$input_file".sorted | sort -u > "$FileName"_distinct_features.txt
