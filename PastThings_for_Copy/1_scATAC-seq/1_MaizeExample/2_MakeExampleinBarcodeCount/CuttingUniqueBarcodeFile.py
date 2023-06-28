#chr1    216     217     CB:Z:CCAATGATCCAGAGAG-1_A619_2  +
import sys, glob
import string

SampleName = sys.argv[1] #A619
## Cell Meta data after clustering :
MetaData="/scratch/sb14489/3.scATAC_flo_Past/5.Socrates/4_CombineAll_AfterD/Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.ALL_CELLs.metadata.txt"
# First column: CB:Z:GTCACCTAGCCTGTAT-1_A619_2
## Will do only in Cluster1
List = []
infile = open(MetaData,"r")
for sLine in infile:
    sList = sLine.strip().split("\t")
    #print(sList[0])
    CellId = sList[0].split("-")[0] #CB:Z:TGTACAGAGACTAGGC-4
    Sample = "_".join(sList[0].split("-")[1:len(sList)])
    Cluster = sList[len(sList)-3]
    #print(Cluster)
    #print(SampleName)
    #print(Cluster)
    #print(Sample)
    #print(sList[0])
    if SampleName == Sample and Cluster == str(1):
        Barcode = CellId.replace("CB:Z:","")
        List.append(sList[0])
#print(List)
infile.close()

DataPath="/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"
Outpath="/scratch/sb14489/3.scATAC/1.MaizeExample/4.Bam_FixingBarcode/"
infile = open(DataPath+SampleName+"_Unique.bed","r")
outfile = open(Outpath+"Ex_"+SampleName+"_Unique.bed","w")

#chr1    333     334     CB:Z:ACCATCCTCAGGTCTA-1_A619    -
#chr1    333     334     CB:Z:TAACGGTCACTGTCAA-1_A619    -

for sLine in infile:
    sList = sLine.strip().split("\t")
    if sList[3] in List:
        outfile.write(sLine)

infile.close()
outfile.close()
