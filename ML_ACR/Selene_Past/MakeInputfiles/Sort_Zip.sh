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


while getopts "b:" opt; do
  case $opt in
    b)
      input_file="$OPTARG"
      FileName="${input_file%.bed}"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

sort -k1V -k2n -k3n "$input_file" > "$input_file".sorted

bgzip -c "$input_file".sorted > "$input_file".sorted.gz

tabix -p bed "$input_file".sorted.gz

cut -f 4 "$input_file".sorted | sort -u > "$FileName"_distinct_features.txt
