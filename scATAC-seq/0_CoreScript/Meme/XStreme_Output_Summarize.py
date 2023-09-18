
OutputXSteme="/scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/IM_OC_dACR_JASPARMotif_Xstream/streme_out/sequences.tsv"
#motif_ID        motif_ALT_ID    motif_P-value   seq_ID  seq_Score       seq_Class       is_holdout?
#1-AGRGARRGAGAGAGA       STREME-1        2.1e-025        chr9:122213355-122213590        20.87   tp      0
#1-AGRGARRGAGAGAGA       STREME-1        2.1e-025        chr7:10333751-10334308  20.87   tp      0
#1-AGRGARRGAGAGAGA       STREME-1        2.1e-025        chr10:4279363-4280396   20.87   tp      0

def Open_Peak_bed(InfileName):
    infile = open(InfileName,"r")
    List = []
    for sLine in infile:
        sList = sLine.strip().split("\t")
        Pos = sList[0]+":"+sList[1]+"-"+sList[2]
        List.append(Pos)
    return(List)

infile = open(OutputXSteme,"r")
infile.readline()

Dic_byMotif={}
for sLine in infile:
    sList = sLine.strip().split("\t")
    if len(sList)>2:
        Dic_byMotif[sList[1]].append(sList[3])
infile.close()

A619Higher = Open_Peak_bed("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC.A619Higher.Bed")
Bif3Higher = Open_Peak_bed("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC.Bif3Higher.Bed")

nA619Higher = 0
nBif3Higher = 0
print(Dic_byMotif.keys())

for i in Dic_byMotif["STREME-1"]:
    if i in A619Higher:
        nA619Higher+=1
    elif i in Bif3Higher:
        nBif3Higher+=1
print(len(Dic_byMotif["STREME-1"]))
print(nA619Higher)
print(nBif3Higher)
