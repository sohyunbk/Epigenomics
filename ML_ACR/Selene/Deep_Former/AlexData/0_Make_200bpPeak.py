infile = open("/scratch/sb14489/8.ML_ACR/1.InputBed/sorted_Seedling_Peaks.bed","r")
#chr1	34427	34927
#chr1   500 1000
#chr1    650 850 +150 -150
## Or make gap as 204
#chr1   648 852 +148 -148
outfile1 = open("/scratch/sb14489/8.ML_ACR/1.InputBed/sorted_Seedling_Peaks_200bp.bed","w")
outfile2 = open("/scratch/sb14489/8.ML_ACR/1.InputBed/sorted_Seedling_Peaks_204bp.bed","w")

for sLine in infile:
    sList = sLine.strip().split("\t")
    nStart = int(sList[1])
    nEnd = int(sList[2])
    outfile1.write(sList[0]+"\t"+str(nStart+150)+"\t"+str(nEnd-150)+"\t"+sList[3]+"\n")
    outfile2.write(sList[0]+"\t"+str(nStart+148)+"\t"+str(nEnd-148)+"\t"+sList[3]+"\n")
infile.close()
outfile1.close()
outfile2.close()
