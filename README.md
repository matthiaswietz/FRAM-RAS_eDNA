## eDNA metabarcoding of microbial communities  

This repo describes the bioinformatic procedure from raw fastq reads to ASV counts, applying cutadapt primer clipping and DADA2 ASV generation. Sequence files were obtained from amplicon-PCRs using 16S and 18S rRNA primers, sequencing using Illuimna MiSeq technology. In DADA2, each sequencing run was processed individually per the programmer's guidelines (https://benjjneb.github.io/dada2/tutorial_1_8.html), and then merged before chimara removal and taxonomic assignment. These scripts describe the procedure on samples obtained by autonomous Remote Access Samplers (RAS) in the FRAM-LTER. Scripts are however generalizable to any ribosomal metabarcoding of eDNA. 

The repo is organized in the DAD2 processing scripts, separeted by 16S (bac_processing) and 18S (euk_processing) rRNA gene sequencing. Each folder contains the separate DAD2 fril fomr different Illunmiane runs (identifed by five-character ID), and a script MergeChimTax.R for merging sequence tables, chimera removal and taxnomoy assignment.

The top-level directory contains files that process the 16S/18S ASV tables further by remoing negative control counts, refomartin taxonomic rankes if appropriate, and connecting with environemtnal data. The script then calculates alpha-diverity indices on ASV tables. Finally , script lists how the full ASV and metadata tables are subsetted for individual studies, e.g. Wietz et al. (https://www.nature.com/articles/s43705-021-00074-4) and Priest et al. (https://www.biorxiv.org/content/10.1101/2022.08.12.503524v2).

