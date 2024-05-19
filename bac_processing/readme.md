DADA2 scripts to process 16S rRNA amplicon sequences; separately per MiSeq run (identified by five-character ID). 

Primer clipping is done with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) via [this script](./cutadapt.sh).

Individual sequence tables are then merged, chimera-checked and taxonomy-asigned using script [MergeChimTax.R](./MergeChimTax.R).

The complete [ASV table](../output/bacSeqtab.txt), [taxonomy table](../output/bacTax.txt) and [ASV sequences](../output/bacUniques.fasta) are then further processed in script [DataLoad.R](../DataLoad.R).
