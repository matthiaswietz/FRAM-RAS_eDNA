## eDNA metabarcoding of microbial communities  

This repo describes the bioinformatic procedure from raw fastq reads to ASV counts, applying cutadapt primer clipping and DADA2 ASV generation. Sequence files were obtained from amplicon-PCRs using 16S and 18S rRNA primers, sequencing using Illuimna MiSeq technology. In DADA2, each sequencing run was processed individually per the programmer's guidelines (https://benjjneb.github.io/dada2/tutorial_1_8.html), and then merged before chimara removal and taxonomic assignment. These scripts describe the procedure on samples obtained by autonomous Remote Access Samplers (RAS) in the FRAM-LTER. Scripts are however generalizable to any ribosomal metabarcoding of eDNA. 

