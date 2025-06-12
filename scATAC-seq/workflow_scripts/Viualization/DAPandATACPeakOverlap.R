library(GenomicRanges)

ATACpeak <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks.bed")
DAPpeak <- read.table("/scratch/sb14489/7.DAPorChIP/DAPseq_WUS/HB67_WUS1_B73v5_Q30_qval5_finalBl/HB67_WUS1_B73v5_Q30_qval5_finalBl.GEM_events.narrowPeak")

atac_gr <- GRanges(seqnames = ATACpeak$V1, 
                   ranges = IRanges(start = ATACpeak$V2, end = ATACpeak$V3))

dap_gr <- GRanges(seqnames = DAPpeak$V1, 
                  ranges = IRanges(start = DAPpeak$V2, end = DAPpeak$V3))
overlaps <- findOverlaps(dap_gr,atac_gr)
overlapping_dap_peaks <- DAPpeak[queryHits(overlaps), ]
