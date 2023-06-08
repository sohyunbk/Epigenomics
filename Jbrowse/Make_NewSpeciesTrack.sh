## Sapelo2
## Fa file and giff should have same chromosome name

conda activate Jbrowse
/home/hj43453/condaenv/test-env/bin/prepare-refseqs.pl --fasta /scratch/hj43453/workDir/02_methylome/10_Poplar/methylpy_index/g717v4_poplar_h1h2_ChrC.fa --out poplar_g717v4
bin/flatfile-to-json.pl --gff /scratch/hj43453/workDir/02_methylome/10_Poplar/methylpy_index/annotation_v4/g717_poplar_genes_h1h2_v4_browserHJ.gff --trackLabel genes --out poplar_g717v4
bin/generate-names.pl -v --out poplar_g717v4


#ssh schmitzlab1@heredity.genetics.uga.edu
#pw: Schmacct5$
#cd /data01/epigenome/JBrowse/
