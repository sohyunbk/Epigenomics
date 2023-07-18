##### This is after RemoveMulitiMap_Deduplication.sh
## From pablo # Sohyun edited for cellranger v2 and do not consider non nuclear configs
#Use : python /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/FindTn5Insertion_fromBam.py  \
#-BAM //scratch/sb14489/3.scATAC/2.Maize_ear/3.SortedBam/"${SampleNameList[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam \
#-exp_name "${SampleNameList[SLURM_ARRAY_TASK_ID]}" -threads 10 \
#-output_file /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${SampleNameList[SLURM_ARRAY_TASK_ID]}"

import pysam
import argparse
import ast
import os
import re
import sys
import numpy as np
import concurrent.futures

def process_read(read):
    try:
        original_tag = read.get_tag("CB")
        exp_name_tag = "-" + exp_name
        new_tag = original_tag.replace("-1", exp_name_tag)
        read.set_tag("CB", new_tag, replace=True)
        return read
    except KeyError:
        return None

def read_bam_file(bam_file, exp_name, outputName, num_threads):
    read_bam_file = pysam.AlignmentFile(bam_file, "rb")
    outfile = pysam.AlignmentFile(outputName+"_BarcodeFixed.sam", "wh", template=read_bam_file.header)

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        processed_reads = executor.map(process_read, read_bam_file)

    for read in processed_reads:
        if read is not None:
            outfile.write(read)

    read_bam_file.close()
    outfile.close()

def makeTn5bed(input_sam_fl, output_dir):
    # Rest of the code for makeTn5bed function here
    ##transfer the
    ##prepare the flag_dic
    flag_dic = {'0':'read_paired',
                '1':'read_properly',
                '2':'read_umapped',
                '3':'mate_unmapped',
                '4':'read_reverse',
                '5':'mate_reverse',
                '6':'first_pair',
                '7':'second_pair',
                '8':'secondary',
                '9':'fail_qc',
                '10':'duplicate',
                '11':'sup_align'
    }

    store_final_line_list = []

    #with open (input_sam_fl, 'r') as ipt:
    opt = open(output_dir+"_Unique.bed","w")
    for eachline in input_sam_fl:
        col = eachline.strip().split() #['A00600:145:HCVY7DSX2:2:2524:25265:24048', '163', 'chr1', '11483', '31', '131M', '=', '11483', '131', 'GTGTACGAGCCTCTGGTCGATGATCAATGGCCACACAACCCCCAATTTTTATGAAAATAGCCATGAGAGACCATTTTCAATAATACTAGAGGCTAAGACCTACAGATTTTTGACCAAGAAATGGTCTCCAC', 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF', 'BC:Z:GCCTCCGT', 'MC:Z:131M', 'MD:Z:103G27', 'PG:Z:MarkDuplicates', #'RG:Z:3_bif3:MissingLibrary:1:HCVY7DSX2:2', 'NM:i:1', 'GP:i:11482', 'MP:i:11613', 'MQ:i:31', 'TQ:Z:FFFFFFF,FFFFFFFFFFF:', 'CR:Z:GAACTTGGTTTAGAAG', 'TR:Z:CTGTCTCTTATACACATCTG', 'AS:i:126', 'XS:i:123', 'QT:Z:FFFFFFFF', 'CY:Z:FFFFFFFFFFFFFFF,', 'CB:Z:GAACTTGGTTTAGAAG-Ex']

        #print(col)
        ##extracct barcode information
        for i in range(9,len(col)):
            if col[i].startswith('CB:Z:'): #CB:Z:AGTTTGGTCCAAACCA-Ex
                bc = col[i]

        ##store flag list information
        flag_list = [] ##[read_reverse,second_pair,duplicate,sup_align]

        ##transfer the flag to binary information
        bin = np.binary_repr(int(col[1]), width=12)
        #print(bin)

        bin_list = list(bin) #['0', '0', '0', '0', '0', '1', '1', '0', '0', '0', '1', '1']
        ##generate flag information
        for i in range(len(bin_list)):
            if bin_list[i] == '1':
                flag_list.append(flag_dic[str(i)])

        #print(flag_list) ['mate_reverse', 'first_pair', 'duplicate', 'sup_align']
        ##get start and end pos of read
        chr = col[2]
        pos1 = col[3]
        cigar = col[5] #131M or 79S72M
        netdif = 1
        #print(cigar)
        ##in order to make act list we need to do transfer the CIGAR string to another format:
        ##dic_list is eg. [{'M': 76}, {'I': 15}, {'M': 57}, {'S': 3}]
        list_CIGAR = re.findall('\d+|\D+', cigar) #['148', 'M']
        #print(list_CIGAR)
        n = int(len(list_CIGAR) / 2)
        dic_list = []
        for i in range(0, int(n)):
            dic = {list_CIGAR[(2 * i + 1)]: int(list_CIGAR[(2 * i)])}
            #print(dic)
            dic_list.append(dic) # it just transfer ['72', 'M'] to {'M': 72}

        cig_dic = {} ####cig_dic stores key and value. key is the 1_M and value is number besides the key in the CIGAR eg {'1_M':50,'2_S':30}
        act_list = [] ##transfer the cigar to [1_M,2_S] or others [1_S,2_M,3_S] stored in the act_list
        ##generate act_list
        item_count = 0
        for eachdic in dic_list:
            item_count += 1
            item_str = str(item_count) + '_' + list(eachdic.keys())[0]
            act_list.append(item_str)
            cig_dic[item_str] = str(eachdic[list(eachdic.keys())[0]])
        #print(act_list)
        for eachact in act_list:
            ##eachact is 1_M or 2_S or others
            ##eachact_list = [1,M] or [2,S]
            ##eachact_list[1] is 'M'
            eachact_list = eachact.split('_')

            ##save the time value to the cig dictionary
            ##time indicates the value in the CIGAR eg. 50M. time is '50'
            time = cig_dic[eachact]

            if eachact_list[1] == 'M' or eachact_list[1] == 'S':
                netdif += int(time)

            elif eachact_list[1] == 'D':
                netdif += int(time)

        #print("End Net difference")
        #print(netdif)
        #print("Updated position")
        pos2 = netdif + int(pos1)

        ##shift
        if flag_list[0] == 'read_reverse':
            end = pos2 - 4
            start = end - 1
            #print (chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '-')
            final_line = chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '-'
            ## Write final line
            opt.write(final_line + '\n')

        else:
            start = int(pos1) + 5
            end = start + 1
            #print (chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '+')
            final_line = chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '+'
            ## Write final line
            opt.write(final_line + '\n')
def get_parser():
    parser = argparse.ArgumentParser(description='Pull out reads aligning to a \
        region from multiple list of BAM files, puts them into a BAM file \
        for later assembly.')
    parser.add_argument('-BAM', "--bam_file", help="Bam file to pull reads from.", required=True, dest='bam_f')
    parser.add_argument('-exp_name', "--experiment_name", help="10x config file to pull scaffold names of nuclear and non-nuclear scaffolds", required=True, dest='exp')
    parser.add_argument('-output_file', "--output", help="Output is Path + prefix.", required=False, dest='o')
    parser.add_argument('-threads', "--num_threads", type=int, help="Number of threads to use for processing (default: 1).", default=1)
    args = vars(parser.parse_args())
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()

    read_bam_file(args.bam_f, args.exp, args.o, args.num_threads)
    makeTn5bed(BarcodedFixedSamFileName, args.o)
