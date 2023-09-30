import pybedtools
from pybedtools import BedTool
import glob
Path = "/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/"
IMOCFile = "IM-OC_FDR.0.01_Bif3Higher.Bed"

IMOC = BedTool(Path+IMOCFile)

all_files = glob.glob(Path+"*.FDR0.05.Bed")
excluded_files = [Path+'IM-OC.FDR0.05.Bed']
filtered_files = [file for file in all_files if file not in excluded_files]

for Other_dACR in filtered_files:
    OtherCellT_dACR = BedTool(Other_dACR)
    Intersect = OtherCellT_dACR.intersect(IMOC)
    IMOC = IMOC - Intersect

IMOC.saveas(Path+IMOCFile.replace(".Bed","_NotOverlapWithOtherCTdACRs.Bed"))

AnnoBed = "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed"
Gene = BedTool(AnnoBed)
AnnIntersect = Gene.intersect(IMOC)
IMOC_NotOverlapTSSorGeneEnd = IMOC-AnnIntersect
IMOC_NotOverlapTSSorGeneEnd.saveas(Path+IMOCFile.replace(".Bed","_NotOverlapWithOtherCTdACRs_NotOverlapTSSorGeneEnd.Bed"))
