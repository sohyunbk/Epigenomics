import os, glob

def CopySummit(InputDir,OutputDir):
  for sDir in os.listdir(InputDir):
      cmd = "cp "+InputDir+sDir+"/*.reproducible_summits.passing_FDR "+OutputDir
      print(cmd)
      os.system(cmd)


CopySummit("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/","/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619")
CopySummit("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3/","/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/Bif3")
