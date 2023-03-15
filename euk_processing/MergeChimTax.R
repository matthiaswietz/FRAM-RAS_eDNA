## Merge all separate runs ##
## Chimera removal, taxonomy with PR2 ##

setwd("/isibhv/projects/FRAMdata/FRAM_MicrObs/WaterCol_RAS_1620/eukaryotes/")
require(dada2)

# Load individual seqtabs 
sq1 <- readRDS(
  "./BYTKM/seqtab_euk_BYTKM.rds")
sq2 <- readRDS(
  "./C947H/seqtab_euk_C947H.rds")
sq3 <- readRDS(
  "./CW5KV/seqtab_euk_CW5KV.rds")
sq4 <- readRDS(
  "./CVD36/seqtab_euk_CVD36.rds")
sq5 <- readRDS(
  "./CG5TF/seqtab_euk_CG5TF.rds")
sq6 <- readRDS(
  "./JK6MD/seqtab_euk_JK6MD.rds")
sq7 <- readRDS(
  "./JR82J/seqtab_euk_JR82J.rds")
sq8 <- readRDS(
  "./K3N66/seqtab_euk_K3N66.rds")
sq9 <- readRDS(
  "./NegCtr/seqtab_euk_NegCtr.rds")

# merge all years/seq-runs
seqtab <- mergeSequenceTables(
  sq1,sq2,sq3,sq4,
  sq5,sq6,sq7,sq8,sq9, repeats="error")

# Removing chimeras 
# 32674 bimeras of 59557 sequences
# allnew: 32025 of 53982 sequences
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", 
  multithread=24, verbose=T)

# stats 
# >90% kept - excellent!
dim(seqtab.nochim)  # 28240 seqs -- 385 samples #allnew: 21957 - 375 
summary(rowSums(seqtab.nochim)/rowSums(seqtab))

# Determine read lengths/size range of amplicons
table(rep(nchar(colnames(seqtab.nochim)), 
          colSums(seqtab.nochim)))

# Remove singletons and junk sequences
# "c" adjusted to size range of amplicons
seqtab.nochim2 <- seqtab.nochim[, nchar(
  colnames(seqtab.nochim)) %in% c(300:430) & 
    colSums(seqtab.nochim) > 1]

# Stats
# 21829 sequences -- allnew: 19720
dim(seqtab.nochim2) 
summary(rowSums(seqtab.nochim2))
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))

################################

## TAXONOMY -- PR2 ##
# different db-versions used
# v12 (as in 1st RAS paper)
# newer v13 and v14

tax_v12 <- assignTaxonomy(
  seqtab.nochim2, 
  "../../tax_db/pr2_version_4.12.0_18S_dada2.fasta.gz", 
  tryRC = TRUE,
  multithread = 24)

tax_v13 <- assignTaxonomy(
  seqtab.nochim2, 
  "../../tax_db/pr2_version_4.13.0_18S_dada2.fasta.gz", 
  tryRC = TRUE,
  multithread = 24)

tax_v14 <- assignTaxonomy(
  seqtab.nochim2, 
  "../../tax_db/pr2_version_4.14.0_SSU_dada2.fasta.gz", 
  tryRC = TRUE,
  multithread = 24)

################################

# Export desired version
# Used v12 for comparability with 16-17 & bloom paper
tax <- tax_v12

# Eukaryota 21829 -- allnew: 19720   
summary(tax)

# Check NA on division level (= phylum) 
table(tax[, 1])   
sum(is.na(tax[, 3]))   #2492 allnew: 1985

# really high... true signals or not?
# left all in for now; influence of phylum-NA removal tested later in R 
# tax.good <- tax[!is.na(tax[, 3]),]
# seqtab.nochim2.good <- seqtab.nochim2[, rownames(tax.good)]
# summary(rowSums(seqtab.nochim2.good))

# Format tables
seqtab.nochim2.print <- t(seqtab.nochim2)
all.equal(rownames(
  seqtab.nochim2.print), rownames(tax)) #TRUE

# Add rownames to seqtab
rownames(seqtab.nochim2.print) <- paste(
  "asv", 1:ncol(seqtab.nochim2), sep = "")

# Add rownames to taxtable
tax.print <- tax
rownames(tax.print) <- rownames(
  seqtab.nochim2.print)

# Export
write.table(
  seqtab.nochim2.print,
  "eukSeqtab.txt", 
  sep="\t", quote=F)
write.table(
  tax.print,
  "eukTax.txt", 
  sep="\t", quote=F)
uniquesToFasta(
  seqtab.nochim2,
  "eukUniques.fasta")

# save taxdata
save.image("RAS_eukTax.Rdata")
