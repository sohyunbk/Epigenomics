## This script is starting from the input file "SampleName_CT.reproducible_summits.passing_FDR"
## SampleName_CT.reproducible_summits.passing_FDR how does it look like?
## chr1	435379	435880	chr1__435379__435880__7.736448127352201
## Last column number is normalized Tn5 Score?
import pybedtools
import glob
import os
import argparse
import sys


def remove_overlapping_peaks(py_bed_tool):
    """TODO: Docstring for remove_overlapping_peaks.
    :returns: TODO
    """
    def ID_correct_interval(arg1):
        """TODO: Docstring for ID_correct_interval.
        :arg1: TODO
        :returns: TODO
        """
        split_to_lists = arg1[3].split(",")
        split_intervals_to_string = [i.split("__") for i in split_to_lists]
        sorted_lists = sorted(
            split_intervals_to_string, key=lambda x: max(x[-1]), reverse=True
        )
        most_signifigant = sorted_lists[0]
        sig_intersect = regen_interval_name(most_signifigant)
        generated_interval = pybedtools.create_interval_from_list(sig_intersect)
        return generated_interval
    final_list = []
    for i in py_bed_tool:
        if "," in i[3]:
            most_sig_interval = ID_correct_interval(i)
            final_list.append(most_sig_interval)
        elif "," not in i[3]:
            final_list.append(i)
    generated_bed_file = pybedtools.BedTool(final_list)
    return generated_bed_file

def regen_interval_name(list_1):
    """Regenerate bed interval with correct score after splitting.
    :arg1: TODO
    :returns: TODO
    """
    generate_string = list_1[0] + "__" +  list_1[1] + "__" +  list_1[2] + "__" +  list_1[3]
    final_list = list_1
    final_list[3] = generate_string
    return(final_list)

def process_files(base_dir):
    # Pattern to match all '.reproducible_summits.passing_FDR' files in subdirectories
    pattern = os.path.join(base_dir, '**/*reproducible_summits.passing_FDR')
    # Recursively search for files matching the pattern
    BedFile = pybedtools.BedTool(glob.glob(pattern, recursive=True)[0])
    #for filename in glob.glob(pattern, recursive=True)[1:]:
    #    print(f"Processing file: {filename}")
    #    BedFile.cat(filename,postmerge=False)
    #    print(len(BedFile))
    AllBed = BedFile.cat(*glob.glob(pattern, recursive=True)[1:],postmerge=False)
    return(AllBed)

def remove_invaid_merge(feature):
    """ Given a bedtool feature, go through and replace the feature name. All
    Features will have the same name.
    """
    if feature.start < feature.stop and feature.start > 0:
        return feature
    elif feature.start < feature.stop and feature.start < 0:
        pass
    elif feature.start > feature.stop:
        pass
    else:
        pass

def extend_fields(feature, n):
    fields = feature.fields[:]
    while len(fields) < n:
        fields.append(".")
    return pybedtools.create_interval_from_list(fields)

def add_final_counter(bed_tool, base_name):
    """TODO: Docstring for add_final_counter.
    :returns: TODO
    """
    final_bed_tool = []
    counter = 1
    for i in bed_tool:
        new = i
        ## Grab the Normalized Score
        grabbed_score = i[3].split("_")[-1]
        generate_new_name = base_name + "_" + str(counter)
        new[3] = generate_new_name
        new[4] = grabbed_score
        final_bed_tool.append(new)
        counter += 1
    return final_bed_tool
def get_parser():
    ##
    parser = argparse.ArgumentParser(
        description="Call Peaks for scATAC data. \
    Requires cluster annnotations, as well as BED file ipput."
    )
    parser.add_argument(
        "-Path1_forCellPeak",
        "--Path1_forCellPeak",
        help="Path1_forCellPeak",
        required=True,
        dest="Path1",
    )
    parser.add_argument(
        "-Path2_forCellPeak",
        "--Path2_forCellPeak",
        help="Path2_forCellPeak",
        required=True,
        dest="Path2",
    )
    parser.add_argument(
        "-OutPutFile", "---OutPutFile", help="Output File. Dir should be made",
         required=False, dest="OutFileName"
    )
    args = vars(parser.parse_args())
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    ## It's actually pretty easy code.
    base_dir1 = args.Path1
    base_dir2 = args.Path2
    OutputFileName = args.OutFileName
    #base_dir1 = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/rel2"
    #base_dir2 = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619"
    #OutputFileName = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed"
    All_bed1 = process_files(base_dir1)
    All_bed2 = process_files(base_dir2)

    ## 1. Combine all the peak files.
    CombineAll_bed = All_bed1.cat(All_bed2,postmerge=False)

    cated_bed_files_cleaned = CombineAll_bed.each(remove_invaid_merge).sort()
    ## Merge peaks from different cell  types select most signifigant

    ## 2. Merge the peaks overlapped but keep the last column as it has cell type name and choose only the peak with the highest Tn5
    merged_cleaned_bed_files = cated_bed_files_cleaned.merge(c=4, o="collapse").each(remove_invaid_merge)

    print(f"Finishing Up Peak Identification")
    generate_final_peak_set = remove_overlapping_peaks(merged_cleaned_bed_files)
    numbered_acr_peaks = pybedtools.BedTool(generate_final_peak_set).sort()
    number_acr_fields = numbered_acr_peaks.field_count() + 1
    numbered_acr_peaks_extended = numbered_acr_peaks.each(
        extend_fields, number_acr_fields
    )
    ## Adding counter to final ACR names
    numbered_acr_peaks_final = add_final_counter(numbered_acr_peaks_extended, "scACR")
    with open(OutputFileName, 'w') as file: file.write(''.join(str(item) for item in numbered_acr_peaks_final))
