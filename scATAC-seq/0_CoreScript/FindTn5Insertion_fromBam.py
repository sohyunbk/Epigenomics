##### This is after RemoveMulitiMap_Deduplication.sh

import pysam
import argparse
import ast
import os
import re
import sys
import numpy as np

## From pablo # Sohyun edited for cellranger v2 and do not consider non nuclear configs
## Usage: python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/4_BarcodeArrange/4-1_FixingBarcodeName.py \
# -BAM "${List[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam -exp_name Ex | samtools view -@ 12 -h - > ../4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_BarcodeFixed.sam

def read_bam_file(bam_file,exp_name,BarcodedFixedSamFileName):
    """
    First alter the CB tags as well as the other tags
    Next - count the tag and add tag to dictionary if not presenat and start at
    one. Add to total column, add to nuclear/nonnuclear column depending on the
    scaffold name and the list given above
    """
    #Read the File
    #save = pysam.set_verbosity(0)
    #read_bam_file = pysam.AlignmentFile(bam_file,"rb", ignore_truncation=True)
    #pysam.set_verbosity(save)
    read_bam_file = pysam.AlignmentFile(bam_file,"rb")

    outfile = pysam.AlignmentFile(BarcodedFixedSamFileName, "wh", template=read_bam_file.header)
    for read in read_bam_file :
        #For each read alter CB tag
        try:
            original_tag = (read.get_tag("CB"))
            #print("here")
            #print(original_tag)
            exp_name_tag = "-" + exp_name
            new_tag = original_tag.replace("-1", exp_name_tag)
            read.set_tag("CB", new_tag, replace=True)
            #set_tag: https://pysam.readthedocs.io/en/latest/api.html?highlight=set_tag#pysam.AlignedSegment.set_tag
            outfile.write(read)
        except KeyError:
            pass

def makeTn5bed(input_sam_fl,output_dir):

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
    opt = open(output_dir,"w")
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
    parser = argparse.ArgumentParser(description='Pull our reads aligning to a\
        region from multiple list of BAM files, puts them into a BAM file\
        for later assembly.')
    parser.add_argument('-BAM', "--bam_file", help="Bam file to \
        pull reads from.", required=True, dest='bam_f')
    parser.add_argument('-exp_name', "--experiment_name", help="10x config file to \
        pull scaffold names of nulcear and non-nuclear scaffolds",
        required=True, dest='exp')

    parser.add_argument('-output_file', "--output", help="Output file to write to. \
        If none given output writes to sOutput file to write to. If none \
        given output writes to sout.", required=False, dest='o')
    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    #Load all Bed files
    BarcodedFixedSamFileName = args.bam_f.replace(".bam","_BarcodeFixed.sam")
    gathered_read_dict = read_bam_file(args.bam_f, args.exp,BarcodedFixedSamFileName)
    makeTn5bed(BarcodedFixedSamFileName, args.o)
