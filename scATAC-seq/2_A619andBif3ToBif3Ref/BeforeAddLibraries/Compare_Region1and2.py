from pybedtools import BedTool
import pybedtools
import argparse

#(1) chr2:3679079..3679413
#(2) chr2:3689915..3690478

MetaData = "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
Meta = open(MetaData,"r")
CellID_IMOC = []
for sLine in Meta:
    sList = sLine.strip().split("\t")
    CellID = sList[0]
    AnnV3 = sList[len(sList)-1]
    if AnnV3 == "IM-OC":
        CellID_IMOC.append(CellID)
Meta.close()

def BringBedFileFromCellID(BedFileName,CellIDlist):
    list =[]
    Bed = open(BedFileName,"r")
    for sLine in Bed:
        sList =sLine.strip().split("\t")
        CellID = sList[3]
        if CellID in CellID_IMOC:
            if sLine.startswith("chr"):
                list.append(sLine)
    Bed.close()
    return list

def FilterChr(list):
    NewList =[]
    for i in list:
        if i.startswith("chr"):
            NewList.append(i)
    return NewList

BifRe2_toRef = BringBedFileFromCellID("/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode_withReadName/3_bif3_2_Unique.bed",CellID_IMOC)
BifRe2_toRef_Filtered = FilterChr(BifRe2_toRef)

BifRe2_toBif3 = BringBedFileFromCellID("/scratch/sb14489/3.scATAC/4.Bif3Ref/4.Bam_FixingBarcode_withReadName/3_bif3_2_Unique.bed",CellID_IMOC)

OriginalZmWUS1Promoter = ["chr2\t3679079\t3679413"]
OriginalZmWUS1Promoter = BedTool(OriginalZmWUS1Promoter)

AddedZmWUS1Promoter = ["chr2\t3689915\t3690478"]
AddedZmWUS1Promoter = BedTool(AddedZmWUS1Promoter)

BifRe2_toRef_Bedtools = BedTool(BifRe2_toRef_Filtered)
BifRe2_toBif3_Bedtools = BedTool(BifRe2_toBif3)

Intersect_ZmWUSPromoter = BifRe2_toRef_Bedtools.intersect(OriginalZmWUS1Promoter,wa=True, wb=True)
Intersect_FakePromoter = BifRe2_toBif3_Bedtools.intersect(AddedZmWUS1Promoter,wa=True, wb=True)

a = BifRe2_toRef_Bedtools.saveas('/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/Ref3_Re2_IMOC_toA619Ref.bed')
d = BifRe2_toBif3_Bedtools.saveas('/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/Ref3_Re2_IMOC_toBif3Ref.bed')

b = Intersect_ZmWUSPromoter.saveas('/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/Bif3_Re2_ToA619Ref_ZmWUS1PromoterPeak.intersect')
c = Intersect_FakePromoter.saveas('/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/Bif3_Re2_ToBif3Ref_Added500bp.intersect')

Intersect_ZmWUSPromoter = open('/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/Bif3_Re2_ToA619Ref_ZmWUS1PromoterPeak.intersect',"r")
Intersect_FakePromoter = open('/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/Bif3_Re2_ToBif3Ref_Added500bp.intersect',"r")
OriginalReads = []
for sLine in Intersect_ZmWUSPromoter:
    OriginalReads.append(str(sLine).split("\t")[5])
k = 0
for sLine in Intersect_FakePromoter:
    if str(sLine).split("\t")[5] in OriginalReads:
        k+=1
        print(str(sLine).split("\t")[5])
#Intersect_Peak_Ann = CommonPeak_RemoveBlackList.intersect(Ann,wa=True, wb=True)

A00257:713:HNVFFDSX2:2:1547:15447:6480
A00257:713:HNVFFDSX2:2:2364:17607:1767
A00257:713:HNVFFDSX2:2:1207:20039:19601
A00257:713:HNVFFDSX2:2:2611:20419:14465
A00257:713:HNVFFDSX2:2:2611:20419:14465
A00257:713:HNVFFDSX2:2:1247:22001:25723
A00257:713:HNVFFDSX2:2:1247:22001:25723
A00257:713:HNVFFDSX2:2:1103:32280:22075
A00257:713:HNVFFDSX2:2:1425:21224:13792
A00257:713:HNVFFDSX2:2:1425:21224:13792
A00257:713:HNVFFDSX2:2:1103:32280:22075
A00257:713:HNVFFDSX2:2:1427:19316:35759
A00257:713:HNVFFDSX2:2:1427:19316:35759
A00257:713:HNVFFDSX2:2:1139:31503:21324
A00257:713:HNVFFDSX2:2:1139:31503:21324
A00257:713:HNVFFDSX2:2:2154:4336:26256
A00257:713:HNVFFDSX2:2:2154:4336:26256
A00257:713:HNVFFDSX2:2:1207:20039:19601
A00257:713:HNVFFDSX2:2:1452:24035:16376
A00257:713:HNVFFDSX2:2:1452:24035:16376
A00257:713:HNVFFDSX2:2:1403:25337:22858
A00257:713:HNVFFDSX2:2:1403:25337:22858
A00257:713:HNVFFDSX2:2:2672:18412:17440
A00257:713:HNVFFDSX2:2:2672:18412:17440
A00257:713:HNVFFDSX2:2:1651:4074:14058
A00257:713:HNVFFDSX2:2:1651:4074:14058
