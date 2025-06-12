import pybedtools

# Define the path to your BED file
CommonPeak = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619_bif3_For_dACR/IM-OC_ComA619Bif3_BLRemove.bed"
A619_PeakFile = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/IM-OC/A619_IM-OC.reproducible_summits.passing_FDR"
Bif3_PeakFile = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3/IM-OC/bif3_IM-OC.reproducible_summits.passing_FDR"

# Load the BED file using pybedtools
def ReadingOnlyChr(File):
    List = []
    infile = open(File,"r")
    for sLine in infile:
        if sLine.startswith("chr"):
            List.append(sLine)
    return(List)
    infile.close()

A619_onlyChr = ReadingOnlyChr(A619_PeakFile)
A619_Peak_bed = pybedtools.BedTool(A619_onlyChr)
CommonPeak_bed = pybedtools.BedTool(CommonPeak)

A619_withCommonPeak = CommonPeak_bed.intersect(A619_Peak_bed)

Bif3_onlyChr = ReadingOnlyChr(Bif3_PeakFile)
Bif3_Peak_bed = pybedtools.BedTool(Bif3_onlyChr)

Bif3_withCommonPeak = CommonPeak_bed.intersect(Bif3_Peak_bed)


unique_to_A619 = A619_withCommonPeak.intersect(Bif3_withCommonPeak, v=True)
