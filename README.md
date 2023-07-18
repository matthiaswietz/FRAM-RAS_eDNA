## eDNA metabarcoding of microbial communities  

This repo describes metabarcoding analysis based on environmental DNA (eDNA), collected using autonomous Remote Access Samplers (RAS) in the [FRAM long-term observatory](https://www.awi.de/en/expedition/observatories/ocean-fram.html) of the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research. Scripts are however generalizable to any ribosomal metabarcoding of eDNA. 

The bioinformatic pipeline comprises all steps from raw fastq reads to ASV counts, including primer clipping and generation of amplicon sequence variants (ASVs). Reads were obtained by PCR with 16S and 18S rRNA primers, followed by Illumina MiSeq sequencing and ASV generation with DADA2. Each Illumina run was processed individually following [the DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial_1_8.html), and then merged before chimera removal and taxonomic assignment. 

All raw sequence files are stored under ENA Umbrella [PRJEB43905](https://www.ebi.ac.uk/ena/browser/view/PRJEB43905).

### Organization of directories and files 

- [bac_processing](./bac_processing): Processing of 16S rRNA reads using individual DADA2 scripts per Illumina run, followed by script [MergeChimTax.R](./bac_processing/MergeChimTax.R) for merging sequence tables, chimera removal and taxonomy assignment [(SILVA database)](https://zenodo.org/record/4587955). 

- [euk_processing](./euk_processing): Processing of 18S rRNA reads using individual DADA2 scripts per Illumina run, followed by script [MergeChimTax.R](./euk_processing/MergeChimTax.R) for merging sequence tables, chimera removal and taxonomy assignment [(PR2 database)](https://github.com/pr2database/pr2database/releases). 

- [output](./output): ASV table, taxonony table, and ASV sequences from both 16S and 18S metabarcoding. Furthermore, ENA accession numbers of all raw fastq files for [16S](./output/ENA_16S_fastq.txt) and [18S](./output/ENA_18S_fastq.txt) amplicons. 

- [metadata](./metadata):  physicochemical measurements and general sample information, needed for detailed analyses as described in the following. Original fastq files and sample information can be matched via columns "sample_title" listed in [sample_info.txt](./metadata/sample_info.txt) and ENA txtfiles in the [output](./output) directory.

This top-level directory contains Rscripts to further process the original data. This includes script [DataLoad.R](./DataLoad.R) to account for negative control counts, refomat taxonomic names if appropriate, and connect with environmental data deposited in [metadata](./metadata). In case samples from several timepoints were pooled, a "mean date" is calculated, and the corresponding environmental parameters averaged as well. 

The script [RarefacDiversity.R](./RarefacDiversity.R) then calculates alpha-diversity indices on ASV tables. Finally, script [DataExport.R](./DataExport.R) subsets the full ASV and metadata tables for individual studies, e.g. describing [microbial dynamics over polar day and night](https://www.nature.com/articles/s43705-021-00074-4) and the [Biological Atlantification of ice-covered waters](https://www.nature.com/articles/s41396-023-01461-6).
