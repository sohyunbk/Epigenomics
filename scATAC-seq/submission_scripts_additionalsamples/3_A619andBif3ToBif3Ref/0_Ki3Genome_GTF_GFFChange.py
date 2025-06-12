#conda activate r_env
## it's python
import bioinfokit
from matplotlib_venn import venn2
from matplotlib import ft2font
from bioinfokit.analys import gff
import os
##~/.conda/env/bin/python3.7
os.chdir("/scratch/sb14489/0.Reference/Maize_Ki3")
gff.gff_to_gtf(file="Zm-Ki3-REFERENCE-NAM-1.0_Zm00029ab.1_OnlyChr.gff3")
