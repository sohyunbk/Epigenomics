import sys, glob
import string

SampleName = sys.argv[1] #A619
## Cell Meta data after clustering :
MetaData="/scratch/sb14489/3.scATAC_flo_Past/5.Socrates/4_CombineAll_AfterD/Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.ALL_CELLs.metadata.txt"
# First column: CB:Z:GTCACCTAGCCTGTAT-1_A619_2
## Will do only in Cluster1
def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))
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
    if SampleName == Sample and Cluster == str(1):
        Barcode = CellId.replace("CB:Z:","")
        ReverseBarcode = rev_compl(Barcode)
        #print("Normal")
        #print(Barcode)
        #print("Reverse")
        #print(ReverseBarcode)
        List.append(Barcode)
        List.append(ReverseBarcode)
infile.close()

##########################################################################################
## R2 file Make Dic
Dic ={}
Switch="OFF"
for sFiles in glob.glob("/scratch/sb14489/3.scATAC/1.MaizeExample/1.Rawdata/OriginalData/"+SampleName+"_S*_R2_*"):
    #print(sFiles)
    infile = open(sFiles,"r")
    outfileName = sFiles.replace(SampleName,"Ex_"+SampleName)
    print(outfileName)
    outfile = open(outfileName,"w")
    j = 0
    for sLine in infile:
        j+=1
        if j%4 == 1:
            Key=sLine
        elif j%4 ==2:
            #print(sLine.strip())
            if sLine.strip() in List:
                #print("Check")
                #print(Key)
                Switch = "On"
                Dic.setdefault(Key,[])
                Dic[Key].append(sLine)
                outfile.write(Key)
                outfile.write(sLine)
            else :
                Switch = "OFF"
        elif j%4 ==3:
            if Switch == "On":
                Dic[Key].append(sLine)
                outfile.write(sLine)
        elif j%4 ==0:
            if Switch == "On":
                Dic[Key].append(sLine)
                outfile.write(sLine)
                Switch ="OFF"
                j =0


    infile.close()
    outfile.close()
#########################################################################################
Switch="OFF"
for sFiles in glob.glob("/scratch/sb14489/3.scATAC/1.MaizeExample/1.Rawdata/OriginalData/"+SampleName+"_S*"):
    if "_R2_" not in sFiles:
        infile = open(sFiles,"r")
        outfile = open(sFiles.replace(SampleName,"Ex_"+SampleName),"w")
        j = 0
        for sLine in infile:
            j+=1
            if j%4 == 1:
                if sLine in Dic:
                    Switch = "On"
                    outfile.write(sLine)
                else :
                    Switch = "OFF"
            elif j%4 ==2:
                if Switch == "On":
                    outfile.write(sLine)
            elif j%4 ==3:
                if Switch == "On":
                    outfile.write(sLine)
            elif j%4 ==0:
                if Switch == "On":
                    outfile.write(sLine)
                    j =0

        infile.close()
        outfile.close()
