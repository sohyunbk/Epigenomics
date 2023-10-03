import pybedtools
from pybedtools import BedTool
Path = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/"
CommonPeakFile = Path+"A619.500bp_peaks.FilteredOrgs.bed"
CommonPeak = BedTool(CommonPeakFile)
CellTypes = ["G2_M","L1","BundleSheath_VascularSchrenchyma","ProcambialMeristem_ProtoXylem_MetaXylem","SPM-base_SM-base","IM-OC","IM_SPM_SM","XylemParenchyma_PithParenchyma","ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma","CalloseRelated","FloralMeristem_SuppressedBract","L1atFloralMeristem","PhloemPrecursor"]
outfile = open("/scratch/sb14489/8.ML_ACR/2.MaizeEar/1.InputBed/MaizeEar_500bpPeak_A619_byCT.bed","w")
for CT in CellTypes:
    FileName = Path+CT+"/A619_"+CT+".reproducible_summits_Onlychr"
    CTSummit= BedTool(FileName)
    Overlap = CommonPeak.intersect(CTSummit,wa=True)
    Overlap = Overlap.merge()
    for i in Overlap:
        outfile.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\t"+CT+"\n")
outfile.close()
