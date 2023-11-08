## 231107 ## Change the Seq one more time.
## Found the reason for the mystery! ZmWUS1 promoter sequences from andrea’s lab was different with Maize v5. Andrea’s lab ZmWUS1 promoter sequences was more matched with sequenced reads.
#Regarding the 445bp promoter sequence of ZmWUS1, it is very similar to the inbred line Ki3, which is one of the NAM founders with a high-quality genome sequence readily available. Probably you can use it as reference.

#Query  1        TTCACAGTAAAGATGCTGGAAAAAAATCATTTGTTGT  37
#                |||||||| ||||||||||||||||||||||||||||
#Sbjct  3436158  TTCACAGTGAAGATGCTGGAAAAAAATCATTTGTTGT  3436194 <- "it startswith "1" but python starts with "0"

## Devide fa file by chromosome

infile = open("/scratch/sb14489/0.Reference/Maize_Ki3/ki3_chromosomes/chr2.fa","r")
#outfile = open("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf_Bif3.fa","w")

Switch="Off"
Chr2=""
#for sLine in infile:
    #if sLine.startswith(">chr2"):
    #    Switch="On"
    #elif sLine.startswith(">chr"):
    #    Switch="Off"
    #elif Switch == "On":
        #Chr2 += sLine.strip()
infile.readline()
for sLine in infile:
    Chr2+=sLine.strip()

infile.close()

print(Chr2[3436157:3436194])

NewSeq="AAAGATATCGTCGTATCCCATGATCTTTCGTGTGTCAACTTCACTTGTCTCTCTCCAAAAGATATCGTCGTATCCCATGATCCTTCCTCCCCTCCAAAAGATATCGTCGTATCaCAGTAAAGCGGCCGAAGCTAGAGCCTTCGTGTGTCAACTTCACTTGTCTCTCTCCAAAAGATATCGTCGTATCCCATGATCCTTCCTCCCCTCCCCTCCCGGCGCCAACCTATATCTCACCATGCACCTAGCACGCAGCTACGCGCGCGCGCGCGCTCTCTCTCTCTCTCTCTCTCTCTGCATGCTAGCTAGCTTTCCTCTAGCCTCTAGCTCCTCGGATATGCACCTATGCACCCCGGCCTCCTTATAAACCCTCCTCAATGCCTCTCCCTTTTCCAAGGCAAGGCCAAGGCGGCAAACCCTTCCCACCGGCCTCCTCCTCCTCCTTGGCGCAGACCGGAGAGATCACAGGAGCTCAGGAAGGCCGGTGTGACCAGCTGCTGAGGTCCTTGGCGCAGACCGGAGAGATCACAGGAGCTCAGGAAGGCCGGTGTGACCAGCTGCTGAGGC" ## 119 + 445bp seq
NewChr2=Chr2[0:3436194]+NewSeq+Chr2[3436194:]

## original fasta file had 80 characters per line

outfile = open("/scratch/sb14489/0.Reference/Maize_Ki3/ki3_chromosomes/chr2_Bif3.fa","w")
outfile.write(">chr2\n")
for i in range(0,len(NewChr2),80):
    if i+80 < len(NewChr2):
        outfile.write(NewChr2[i:i+80]+"\n")
    else:
        outfile.write(NewChr2[i:]+"\n")

outfile.close()
#end : TTTAGGGTTTAGGGT


import os
os.system("cd /scratch/sb14489/0.Reference/Maize_Ki3/ki3_chromosomes/")
os.system("cat chr1.fa chr2_Bif3.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa > ../Zm-Ki3-REFERENCE-NAM-1.0_OnlyChr_Bif3.fa")
