## This script needs cause the output from "De_novo_marker_byCluster" cant be used for markger gene list.
## Why? They have redundant markergenes by cell cluster.
# EX;
#chr	start	end	geneID	name	type
#chr2	150737636	150741196	Zm00001eb093200	Zm00001eb093200_pval_0_FloralMeristem_SuppressedBract	FloralMeristem_SuppressedBract
#chr5	24519843	24527390	Zm00001eb220660	Zm00001eb220660_pval_0_FloralMeristem_SuppressedBract	FloralMeristem_SuppressedBract
#chr4	192613899	192616536	Zm00001eb196470	Zm00001eb196470_pval_0_FloralMeristem_SuppressedBract	FloralMeristem_SuppressedBract
#chr4	192613796	192616522	Zm00001eb196460	Zm00001eb196460_pval_0_FloralMeristem_SuppressedBract	FloralMeristem_SuppressedBract
#chr5	221674642	221681091	Zm00001eb257620	Zm00001eb257620_pval_0_FloralMeristem_SuppressedBract	FloralMeristem_SuppressedBract
#chr3	83813927	83840797	Zm00001eb132450	Zm00001eb132450_pval_0_G2_M	G2_M
#chr1	260696094	260705322	Zm00001eb257620	Zm00001eb257620_pval_0_G2_M	G2_M

#!/bin/bash
#!/bin/bash

# Check if the correct number of arguments was provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

# Assign arguments to variables
input_file="$1"
output_file="$2"

# Process the file and sort by the last column
awk '
BEGIN { FS=OFS="\t" }
NR == 1 { print; next } # Print the header and skip to the next record
{
    # Change the 5th column to match the 4th column
    $5 = $4;

    # Create a key based on the fourth column
    key = $4;

    # Concatenate types with a comma for matching keys
    if (key in data) {
        split(data[key], line, FS);
        line[6] = line[6] "," $6;
        data[key] = line[1] FS line[2] FS line[3] FS line[4] FS line[5] FS line[6];
    } else {
        data[key] = $0;
    }
}
END {
    # Print the merged rows
    for (key in data) {
        print data[key];
    }
}' "$input_file" | sort -t$'\t' -k6,6 > "$output_file"

echo "Processed file saved as $output_file"
