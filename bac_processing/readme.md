Directory containing individual DADA2 scripts per Illumina sequencing run (identifed by five-character ID) of 16S rRNA amplicons from eDNA. 

Individual sequence tables are then merged, chimera-checked and taxonomy-asigned using script [MergeChimTax.R](./MergeChimTax.R).

The complete [ASV table](../output/bacSeqtab.txt), [taxonomy table](../output/bacTax.txt) and [ASV sequences](../output/bacUniques.fasta) are then further processed in script [1_DataLoad.R](../1_DataLoad.R).
