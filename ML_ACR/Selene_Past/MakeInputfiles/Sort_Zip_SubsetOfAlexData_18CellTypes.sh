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

output_file="Seedling_18Celltypes.500.bed"

# Check if output file exists and remove it
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# List of specified filenames
#abaxial_bundle_sheath
#adaxial_leaf_primordia
#bundle_sheath
#cortex
#dividing_leaf_primordia
#ground_meristem
#guard_mother_cell
#hypodermal_sclerenchyma
#L1_leaf_primordia_boundary
#leaf_primordia
#mesophyll
#mesophyll_precursors
#phloem_SE_procambial_precursors
#pith_parenchyma
#procambial_meristem
#protodermal_cell
#protophloem_SE
#xylem_parenchyma

files=(
    "abaxial_bundle_sheath.accessible_ACRs.bed.500"
    "adaxial_leaf_primordia.accessible_ACRs.bed.500"
    "bundle_sheath.accessible_ACRs.bed.500"
    "cortex.accessible_ACRs.bed.500"
    "dividing_leaf_primordia.accessible_ACRs.bed.500"
    "ground_meristem.accessible_ACRs.bed.500"
    "guard_mother_cell.accessible_ACRs.bed.500"
    "hypodermal_sclerenchyma.accessible_ACRs.bed.500"
    "L1_leaf_primordia_boundary.accessible_ACRs.bed.500"
    "leaf_primordia.accessible_ACRs.bed.500"
    "mesophyll.accessible_ACRs.bed.500"
    "mesophyll_precursors.accessible_ACRs.bed.500"
    "phloem_SE_procambial_precursors.accessible_ACRs.bed.500"
    "pith_parenchyma.accessible_ACRs.bed.500"
    "procambial_meristem.accessible_ACRs.bed.500"
    "protodermal_cell.accessible_ACRs.bed.500"
    "protophloem_SE.accessible_ACRs.bed.500"
    "xylem_parenchyma.accessible_ACRs.bed.500"
)

# For each specified file
for file in "${files[@]}"; do
    # Extract the name before the first '.' as label
    label=$(echo $file | cut -f 1 -d '.')

    # Append the label as the fourth column and write to the output file
    awk -v label="$label" '{print $1 "\t" $2 "\t" $3 "\t" label}' $file >> $output_file
done

input_file="Seedling_18Celltypes.500.bed"
FileName=Seedling_18Celltypes
conda activate /home/sb14489/miniconda3/envs/pytorch

sort -k1V -k2n -k3n "$input_file" > "$input_file".sorted

bgzip -c "$input_file".sorted > "$input_file".sorted.gz

tabix -p bed "$input_file".sorted.gz

cut -f 4 "$input_file".sorted | sort -u > "$FileName"_distinct_features.txt
