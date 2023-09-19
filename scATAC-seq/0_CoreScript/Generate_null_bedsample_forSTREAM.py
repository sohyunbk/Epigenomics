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

def compute_gc_ratio_OneSequence(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_length = len(sequence)
    gc_ratio = gc_count / total_length
    return gc_ratio

def get_parser():
    parser = argparse.ArgumentParser(description='Finds peaks shared between \
        replicate peak calls, as well as unqiue peaks to each replicate and \
        outputs said peaks. ')
    parser.add_argument('-bed','--bed_file', help='Bed File Target',\
        required=True, dest='bed'),
    parser.add_argument('-genome','--genome_file', help='Genome File to use', \
        required=True, dest='gn'),
    parser.add_argument('-genome_index','--genome_index', help='Genome File to use Do not have MtPt', \
        required=True, dest='gni'),
    parser.add_argument('-AllPeakForControl','--AllPeakForControl', help='Genome File to use Do not have MtPt', \
        required=True, dest='AllPeakForControl'),
    parser.add_argument('-o','--output_name', help='output',
        required=False, dest='o')
    args = vars(parser.parse_args())
    return parser

if __name__ == "__main__":
    #python scripts/gen_null_bed_sample.py
    # -bed {input.classified_acrs} -genome {params.fasta_file}
    # -genome_index {params.fai} -o {params.control_bed    _output_base}
    args = get_parser().parse_args()
    bed = args.bed
    Genome = args.gn
    Genome_fai = args.gni
    all_peak_minus_for_control=args.AllPeakForControl
    ControlFile= args.o
    #Read Files

    #bed = "/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC.FDR0.05.Bed"
    #Genome = "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa"
    #Genome_fai="/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai"
    #all_peak_minus_for_control="/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/Cat_A619_Bif3_500bpPeak.bed"
    #ControlFile="/scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/Control_for_IM_OC.fa"

    substract_peak_bed=all_peak_minus_for_control.replace(".bed",".Substract.bed")
    bed_file = pybedtools.BedTool(bed)

    ### Target Sequence length distribution #####
    all_length= get_interval_lengths(bed_file)
    max_length = max(bed_file, key=interval_length)
    max_length_value = len(max_length)

    windows = pybedtools.BedTool().window_maker(g=Genome_fai, w=max_length_value)
    all_peak_bed = pybedtools.BedTool(all_peak_minus_for_control)
    ControlRegion = windows.subtract(all_peak_bed)
    ControlRegion.saveas(substract_peak_bed)

    ## Target Sequence GC content distribution #####
    substract_peak_bed_GCContent = substract_peak_bed.replace(".bed",".GCRatio")
    if os.path.exists(substract_peak_bed_GCContent):
        infile = open(substract_peak_bed_GCContent,"r")
        mean_cell_type_gc_score = float(infile.readline().strip())
    else:
        print("The file does not exist.")
        Gc_content = capture_gc_content(bed_file,Genome) ## It takes long like 30 min...
        #print(GC_content[1]):
        #chr1	105053790	105055096	0.420368	0.579632	269	398	359	280	0	1306
        take_len = bed_file.field_count()
        mean_cell_type_gc_score = calcualte_mean_GC_content(Gc_content, take_len)
        OutGC = open(substract_peak_bed_GCContent,"w")
        OutGC.write(str(mean_cell_type_gc_score))
        OutGC.close()

    lower_bound = round(mean_cell_type_gc_score * 0.90, 3)
    upper_bound = round(mean_cell_type_gc_score * 1.10, 3)

    ## RandomSampling
    random.seed(42); RandomN = [random.randint(0, len(all_length)*100) for _ in range(len(ControlRegion))]
    Fulfill = 0
    outfile = open(ControlFile,"w")
    for nRandom in RandomN:
        if Fulfill < len(all_length):
            TargetLegnth = len(bed_file[Fulfill])
            ControlRegionPos = ControlRegion[nRandom]
            ControlRegionPos.end = ControlRegionPos.start + TargetLegnth
            ControlRegionPos_Re = pybedtools.BedTool(str(ControlRegionPos), from_string=True)
            Control_seq = ControlRegionPos_Re.sequence(fi=Genome)
            with open(Control_seq.seqfn, 'r') as f:
                f.readline()
                ControlSequence = f.readline().strip()
            Gc_content_Control = compute_gc_ratio_OneSequence(ControlSequence)
            if lower_bound <  Gc_content_Control <upper_bound :
                outfile.write(">"+ControlRegionPos[0]+":"+ControlRegionPos[1]+"-"+ControlRegionPos[2]+"\n")
                outfile.write(ControlSequence+"\n")
                Fulfill+=1

    if (Fulfill+1) < len(all_length):
        print("Sequence is not enough")

    outfile.close()
