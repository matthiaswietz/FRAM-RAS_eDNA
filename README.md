## eDNA metabarcoding of microbial communities  

This repo describes the bioinformatic pipeline from raw fastq reads to ASV counts, including primer clipping and ASV generation. Fastq files were obtained by PCR of eDNA with 16S and 18S rRNA primers, followed by Illumina MiSeq sequencing. In DADA2, each sequencing run was processed individually per the developer's guidelines (https://benjjneb.github.io/dada2/tutorial_1_8.html), and then merged before chimera removal and taxonomic assignment. 

This repo describes processing of eDNA samples from autonomous Remote Access Samplers (RAS) in the FRAM-LTER (https://www.awi.de/en/expedition/observatories/ocean-fram.html). Scripts are however generalizable to any ribosomal metabarcoding of eDNA. 

The repo is organized into DADA2 processing scripts, separated by [16S](./bac_processing) and [18S](./euk_processing) rRNA gene sequencing. Each folder contains the individual DADA2 pipelines per Illumina run, and a script for merging sequence tables, chimera removal and taxonomy assignment. See the [readme](./bac_processing/readme.md) for detailed information.

This top-level directory contains Rscripts that further process the raw 16S/18S ASV tables. This includes script [DataLoad.R](./DataLoad.R) to account for negative control counts, refomat taxonomic names if appropriate, and connect with environmental data deposited in [metadata](./metadata). In case samples from several timepoints were pooled, a "mean date" is calculated, and the corresponding environmental parameters averaged as well. 

The script [RarefacDiversity.R](./RarefacDiversity.R) then calculates alpha-diversity indices on ASV tables. Finally, script [DataExport.R](./DataExport.R) subsets the full ASV and metadata tables for individual studies, e.g. Wietz et al. (https://www.nature.com/articles/s43705-021-00074-4) and Priest et al. (https://www.biorxiv.org/content/10.1101/2022.08.12.503524v2).

The tables XX and XX contain the formatted, cleaned and averaged (if applicable) data. This file lists the ENA accession numbers of all raw fastq files. 

