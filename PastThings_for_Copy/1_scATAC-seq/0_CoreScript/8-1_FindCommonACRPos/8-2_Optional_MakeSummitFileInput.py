import glob, math
#Current file:Mt	9176	9677	Mt__9176__9677__62.73360575720627
#Edited file: Mt	36422	36423	A619_L1atDeterminateLaterOrgan.pool.macs_peak_10	1193.47

for sFiles in glob.glob("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/*.passing_FDR"):
    print(sFiles)
    infile = open(sFiles,"r")
    outfile = open(sFiles+"_Edited","w")
    for sLine in infile:
        sList = sLine.strip().split("\t")
        sStart = (int(sList[1])+int(sList[2]))/2-0.5
        sEnd = (int(sList[1])+int(sList[2]))/2+0.5
        sScore=sList[3].split("_")[len(sList[3].split("_"))-1]
        NewLine = sList[0]+"\t"+str(math.trunc(sStart))+"\t"+str(math.trunc(sEnd))+"\t"+sList[3]+"\t"+sScore+"\n"
        outfile.write(NewLine)

    infile.close()
    outfile.close()
