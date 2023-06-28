import os, glob

BringFilesPath="/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/tracks/Peaks/"
## GitHubInfo
GithubDir="peaks"
# peaks,ComA619Bif3.unique500bpPeaks,ComA619Bif3.unique500bpPeaks_BeforeFilteringFDR,ACR_Peak,,,,,,,,,,
for sFiles in glob.glob(BringFilesPath+"*"):
    FileNamewithE = os.path.basename(sFiles)
    FileName = FileNamewithE.split(".")[0]
    #SampleName = FileName.split("_")[0]
    CellName = "_".join(FileName.split("_")[1:]).split(".")[0]
    cmd = "%s,%s,%s,scATACMutantsV3/%s,,,,,,,,,,"%(GithubDir,FileNamewithE,FileNamewithE,CellName)
    print(cmd)
