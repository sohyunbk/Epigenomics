import os, glob

BringFilesPath="/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3_MergePeakbyCelltypes_Method1/tracks/"
#BringFilesPath="/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn//tracks/"

## GitHubInfo
GithubDir="peaks"
# peaks,ComA619Bif3.unique500bpPeaks,ComA619Bif3.unique500bpPeaks_BeforeFilteringFDR,ACR_Peak,,,,,,,,,,
for sFiles in glob.glob(BringFilesPath+"*"):
    FileNamewithE = os.path.basename(sFiles)
    FileName = FileNamewithE.split(".")[0]
    #SampleName = FileName.split("_")[0]
    CellName = FileName.split(".")[0]
    cmd = "%s,%s,%s,scATACMutantsV3/%s,,,,,,,,,,"%(GithubDir,FileNamewithE,FileNamewithE,CellName)
    print(cmd)
