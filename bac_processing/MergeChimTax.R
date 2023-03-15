## Merge all separate runs ##
## Chimera removal, taxonomy with Silva138

setwd("/isibhv/projects/FRAMdata/FRAM_MicrObs/WaterCol_RAS_1620/bacteria/")
require(dada2)


#################################

## Merge duplicate Miseq Runs
## BLWBG and BLVR5 (96 samples, annual cycle 2016-17)

temp1 <- readRDS(
  "./BLWBG/seqtab_bac_BLWBG.rds")
temp2 <- readRDS(
  "./BLVR5/seqtab_bac_BLVR5.rds")

sq1 <- mergeSequenceTables(
  temp1, temp2, repeats="sum")

########################################

# Load individual seqtabs from other years
sq2 <- readRDS(
  "./J34R3/seqtab_bac_J34R3.rds")
sq3 <- readRDS(
  "./J3672/seqtab_bac_J3672.rds")
sq4 <- readRDS(
  "./J3DYP/seqtab_bac_J3DYP.rds")
sq5 <- readRDS(
  "./JB7R9/seqtab_bac_JB7R9.rds")
sq6 <- readRDS(
  "./JBBW2/seqtab_bac_JBBW2.rds")
sq7 <- readRDS(
  "./JK6MD/seqtab_bac_JK6MD.rds")
sq8 <- readRDS(
  "./JR82J/seqtab_bac_JR82J.rds")
sq9 <- readRDS(
  "./K3N66/seqtab_bac_K3N66.rds")
sq10 <- readRDS(
  "./NegCtr/seqtab_bac_NegCtr.rds")

# Merge all years/seq-runs
seqtab <- mergeSequenceTables(
  sq1,sq2,sq3,sq4,sq5,
  sq6,sq7,sq8,sq9,sq10, repeats="error")

# Remove chimeras 
# org/aphros: 15621 bimeras out of 35418
# new BLV-MWG: 15411 bimeras out of 32682 
# all new: 15071 bimeras out of 31890 
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method = "consensus", 
  multithread=24, verbose=T)

# 384 samples -- 19797 sequences (b4 357 -- 18599)
# new BLV-MWG:  384 17271
# all new: 384 16819
# >90% kept - good!
dim(seqtab.nochim)  
summary(rowSums(seqtab.nochim)/rowSums(seqtab))

# Determine amplicon length/size range 
table(rep(nchar(
  colnames(seqtab.nochim)), 
  colSums(seqtab.nochim)))

# Remove singletons and junk sequences
# "c" adjusted to size range of amplicons
seqtab.nochim2 <- seqtab.nochim[, nchar(
  colnames(seqtab.nochim)) %in% c(354:415) & 
    colSums(seqtab.nochim) > 1]

# Stats
dim(seqtab.nochim2) # 16805 sequences // #15939 new // #15630 allnew
summary(rowSums(seqtab.nochim2))
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))

###################################################################################

## TAXONOMY -- Silva v138 ##

tax <- assignTaxonomy(
  seqtab.nochim2, 
  "../../tax_db/silva_nr_v138_train_set.fa.gz", 
  tryRC = TRUE,
  multithread = 24)

# Bact 16243 -- Arch 436 -- Euk 81 -- NA 45
# new: Bact 15385 -- Arch 433 -- Euk 74 -- NA 47
# allnew/final: Bact 15108 -- Arch 439 -- Euk 47 -- NA 36
summary(tax)

# Remove NA on phylum level
table(tax[, 1])   
sum(is.na(tax[, 2])) #395
tax.good <- tax[!is.na(tax[, 2]),]
seqtab.nochim2.good <- seqtab.nochim2[
  , rownames(tax.good)]
summary(rowSums(seqtab.nochim2.good))

# Format tables
seqtab.nochim2.print <- t(seqtab.nochim2.good)
tax.print <- tax.good
all.equal(rownames(
  seqtab.nochim2.print), 
  rownames(tax.print)) #TRUE
rownames(seqtab.nochim2.print) <- paste(
  "asv", 1:ncol(seqtab.nochim2.good), sep = "")
rownames(tax.print) <- rownames(seqtab.nochim2.print)

# Export
write.table(
  seqtab.nochim2.print,"bacSeqtab.txt", 
  sep="\t", quote=F)
write.table(
  tax.print,"bacTax.txt", 
  sep="\t", quote=F)
uniquesToFasta(
  seqtab.nochim2.good,
  "bacUniques.fasta")

# save taxdata
save.image("RAS_bacTax.Rdata")
