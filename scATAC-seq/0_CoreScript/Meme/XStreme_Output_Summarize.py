
OutputXSteme="/scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/IM_OC_dACR_JASPARMotif_Xstream/streme_out/sequences.tsv"
#motif_ID        motif_ALT_ID    motif_P-value   seq_ID  seq_Score       seq_Class       is_holdout?
#1-AGRGARRGAGAGAGA       STREME-1        2.1e-025        chr9:122213355-122213590        20.87   tp      0
#1-AGRGARRGAGAGAGA       STREME-1        2.1e-025        chr7:10333751-10334308  20.87   tp      0
#1-AGRGARRGAGAGAGA       STREME-1        2.1e-025        chr10:4279363-4280396   20.87   tp      0


infile = open(OutputXSteme,"r")
infile.readline()

Dic_byMotif={}
for sLine in infile:
