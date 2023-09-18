import pybedtools
import argparse
import os
import itertools
import copy
import statistics
import csv
from threading import Thread
import random
from Bio import SeqIO

def read_in_bed_file(bed_file):
    """
    Reads in a bed file using the pybedtools standard appraoch.
    """
    def intersecting_feature(feature, index):
        """TODO: Docstring for intersecting_feature(feature, L.
        :returns: TODO
        """
        return(feature[index]) != '.'

    try:
        read_file = pybedtools.BedTool(bed_file)
    except:
        print("File DOES NOT EXISTS")
        exit(1)

    return read_file


def capture_gc_content(bed_file, genome_file):
    """Calls GC content using pybedtools. Requires the actual .fa genome.
    Returns the output from that
    :arg1: TODO
    :returns: TODO
    """
    seq_profiled = bed_file.nucleotide_content(fi=genome_file).saveas()
    return(seq_profiled)


def calcualte_mean_GC_content(bed_file_nuc, start_index):
    GC_cont_all = [float(z[start_index + 1]) for z in bed_file_nuc]
    mean_gc_cont = statistics.mean(GC_cont_all)
    return(mean_gc_cont)


def calcualte_mean_counts(bed_file_nuc):
    """
    Calculates the mean count of intersections. The last integer here referes
    to the last column in the intersect command where this is called
    """

    mean_intersection = [int(z[-1]) for z in bed_file_nuc]
    calc_mean = statistics.mean(mean_intersection)
    return(calc_mean)


def calcualte_prop(bed_file_nuc):
    """
    Calculates the mean count of intersections. The last integer here referes
    to the last column in the intersect command where this is called
    """

    total_num = bed_file_nuc.count()

    final_num = bed_file_nuc.filter(lambda x: int(x[-1]) >= 1)
    count_final_num = final_num.count()

    proportion = count_final_num/total_num
    return(proportion)





def generate_matching_GC(ControlRegion,all_length, Genome, mean_cell_type_gc_score):
    print(f"Attempt {current_attempt} out of {max_attempts}")
    if current_attempt > max_attempts:
        raise ValueError("Exceeded maximum attempts to sample matching GC content")
    #SHuffle these regions
    grab_equal_features = acr_bed.random_subset(count)
    #Count the number of fields present after this region survives (for use
    #later when calclating GC()
    surv_index = grab_equal_features.field_count()
    #Calc GC content
    rand_GC_content = capture_gc_content(grab_equal_features, genome_file)
    sample_mean_gc_score = calcualte_mean_GC_content(rand_GC_content,  surv_index)
    print(GC_mean)
    print(sample_mean_gc_score)
    lower_bound = round(GC_mean * 0.90,3)
    upper_bound = round(GC_mean * 1.1,3)
    print(lower_bound, upper_bound)
    if lower_bound <= float(sample_mean_gc_score) <= upper_bound:
        return grab_equal_features
    else:
        return generate_matching_GC(acr_bed, count, genome_file, GC_mean, max_attempts, current_attempt + 1)




def keep_nfeatures(feature, n):
    """Given a list of features, keep only-
    n number of fields

    """
    new_feature = feature[:int(n)]
    feature = new_feature

    return feature




def get_parser():
    parser = argparse.ArgumentParser(description='Finds peaks shared between \
        replicate peak calls, as well as unqiue peaks to each replicate and \
        outputs said peaks. ')
    parser.add_argument('-bed','--bed_file', help='Bed File to Mimic',\
        required=True, dest='bed'),
    #parser.add_argument('-TFs','--TF_file', help='TF file in bed',\
    #    required=True, dest='TF_b'),
    #parser.add_argument('-exclusion_beds','--excl_files', help='Bed File to Mimic',\
    #    nargs="+", required=True, dest='exl'),
    parser.add_argument('-genome','--genome_file', help='Genome File to use', \
        required=True, dest='gn'),
    parser.add_argument('-genome_index','--genome_index', help='Genome File to use', \
        required=True, dest='gni'),
    parser.add_argument('-o','--output_name', help='output',
        required=False, dest='o')
    args = vars(parser.parse_args())
    return parser

def interval_length(interval):
    return len(interval)

def get_interval_lengths(bed_file):
    """
    Get a list of lengths for all intervals in a BedTool object.
    Parameters:
    - bed_file: A BedTool object containing intervals.
    Returns:
    - lengths: A list of lengths for all intervals.
    """
    lengths = []
    for interval in bed_file:
        interval_length = len(interval)
        lengths.append(interval_length)
    return lengths

def calculate_gc_ratio(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    total_bases = len(sequence)
    return gc_count / total_bases

if __name__ == "__main__":
    #python scripts/gen_null_bed_sample.py
    # -bed {input.classified_acrs} -genome {params.fasta_file}
    # -genome_index {params.fai} -o {params.control_bed    _output_base}
    args = get_parser().parse_args()

    #
    #Read Files
    bed = "/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC.FDR0.05.Bed"
    Genome = "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa"
    Genome_fai="/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai"
    all_peak_minus_for_control="/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/Cat_A619_Bif3_500bpPeak.bed"
    substract_peak_bed=all_peak_minus_for_control.replace(".bed",".Substract.bed")
    bed_file = pybedtools.BedTool(bed)
    all_length= get_interval_lengths(bed_file)
    max_length = max(bed_file, key=interval_length)
    max_length_value = len(max_length)


    # Create a BedTool from the FAI file to generate windows
    windows = pybedtools.BedTool().window_maker(g=Genome_fai, w=max_length_value)
    all_peak_bed = pybedtools.BedTool(all_peak_minus_for_control)
    ControlRegion = windows.subtract(all_peak_bed)

    # Save the result to the output file
    ControlRegion.saveas(substract_peak_bed)
    ## Start from here!

    Gc_content = capture_gc_content(bed_file,Genome)
    #print(GC_content[1]):
    #chr1	105053790	105055096	0.420368	0.579632	269	398	359	280	0	1306
    take_len = bed_file.field_count()
    mean_cell_type_gc_score = calcualte_mean_GC_content(Gc_content, take_len)



    ##### Sampling the regions!!
    lower_bound = round(mean_cell_type_gc_score * 0.90, 3)
    upper_bound = round(mean_cell_type_gc_score * 1.10, 3)
    selected_regions = []
    desired_num_regions = bed_file.count()

    # Create a list of regions that don't meet GC ratio criteria


    # Iterate through ControlRegion
