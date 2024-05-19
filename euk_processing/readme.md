Directory containing DADA2 scripts for each Illumina sequencing run (identified by five-character ID) to process 18S rRNA amplicon sequences. 

Primer clipping is done using script [cutadapt.sh](./cutadapt.sh).

Individual sequence tables are then merged, chimera-checked and taxonomy-asigned using script [MergeChimTax.R](./MergeChimTax.R).

The complete [ASV table](../output/eukSeqtab.txt), [taxonomy table](../output/eukTax.txt) and [ASV sequences](../output/eukUniques.fasta) are then further processed in script [DataLoad.R](../DataLoad.R).
