############################################################################################
###  FRAM-RAS -- AMPLICON TIME-SERIES  ###
############################################################################################

# Formatting of dada2 ASVs
# Formatting taxonomy
# Contaminant check
# Alpha-Diversity of cleaned ASV table
# Basic metadata: mean date (if bags pooled)
# Also CTD-depth added (needed for Gfbio)
# Nutrients, ice-cover, chlorophyll

# Master-data for all subsequent studies
# Individual papers will subset this table as needed

####################################

# Set working directory
setwd("/Users/mwietz/ownCloud - mwietz@owncloud.mpi-bremen.de/AWI_MPI/FRAM/RAS/ampliconTimeseries")

library(gtools)
library(dplyr)
library(tidyr)
library(solartime)
#load("_RAS_DataLoad.Rdata")


############################################################################################
###  METADATA ###
############################################################################################

# Load metadata
# Remove neg control columns
# Remove outlier (>200m depth)
ENV <- read.csv(
  "./metadata/sample_info.txt", 
  h=T, sep="\t", stringsAsFactors=F, skipNul=T) %>%
  filter(mooring!="NK") %>%
  filter(!RAS_id %in% c(
    "03_2017_F4_1","04_2017_EGC_2")) %>%
  na_if("") %>% dplyr::select(c(
    "RAS_id","lat","locus_tag",
    "date1","date2","date3","date4",
    "mooring_full","mooring")) %>%
  reshape2::melt(id.vars=c(
    "locus_tag","mooring_full",
    "RAS_id","lat","mooring")) %>%
  dplyr::rename(date=value) 

# Calculate daylight
ENV$daylight <- computeDayLength(
  as.Date(ENV$date, format="%Y-%m-%d"), ENV$lat)

# calculate mean date
# back-convert to char for easier merginn
# will be reassigned to date eventually
meanDate <- read.csv(
  "./metadata/sample_info.txt", h=T, sep="\t", 
    stringsAsFactors=F, skipNul=T) %>%
  filter(mooring!="NK") %>%
  mutate_at(vars(date1, date2, date3, date4),
    as.Date, format="%Y-%m-%d") %>%
  na_if("") %>% rowwise %>%
  mutate(date = mean.Date(c(
    date1, date2, date3, date4), na.rm=T)) %>%
  mutate(date = as.character(date)) %>%
 # mutate(mooring_full = as.character(mooring_full))
  as.data.frame()


############################################################################################
###  IMPORT BACTERIAL / 16S DATA ###
############################################################################################

# ASV table contains sequences from 250m depth
# will be filtered out subsequently (other study)

ASV.bac <- read.table(
  "./bac_output/bacSeqtab.txt",
  h = T, sep = "\t",
  check.names=F)

TAX.bac <- read.table(
  "./bac_output/bacTax.txt",
  h = T, sep = "\t")

####################################

### CONTAMINANT CHECK

# Define samples PCR-amplified with 25/30/35 cycles
# info is contained in meanDate DF
NK25 <- meanDate %>% filter(
  cycles=="25" & locus_tag=="16S") %>% pull(clip_id)
NK30 <- meanDate %>% filter(
  cycles=="30" & locus_tag=="16S") %>% pull(clip_id)
NK35 <- meanDate %>% filter(
  cycles=="35" & locus_tag=="16S") %>% pull(clip_id)

# Define corresponding negative controls
ASV.bac$NK25 <- rowMeans(ASV.bac[,c(
  "16S_NK_1_clip","16S_NK_1_clip",
  "NK_S39_clip","NK_S49_clip",
  "neg-control-bact-25c_clip",
  "negativ_control_clip",
  "PS126_NK_25c_bact_clip")])
ASV.bac$NK30 <- rowMeans(ASV.bac[,c(
  "neg-control-bact-30c_clip",
  "PS126_NK_30c_bact_clip")])
ASV.bac$NK35 <- rowMeans(ASV.bac[,c(
  "neg-control-bact-35c_clip",
  "PS126_NK_35c_bact_clip")])

# round
ASV.bac$NK25 <- round(ASV.bac$NK25, 0)
ASV.bac$NK30 <- round(ASV.bac$NK30, 0)
ASV.bac$NK35 <- round(ASV.bac$NK35, 0)

# Subtract negative counts 
asv1 = ASV.bac[NK25] - ASV.bac$NK25
asv2 = ASV.bac[NK30] - ASV.bac$NK30
asv3 = ASV.bac[NK35] - ASV.bac$NK35

# Join everything again
ASV.bac <- cbind(asv1, asv2, asv3)

# Set negative values to zero
ASV.bac[ASV.bac < 0] <- 0

# Remove mitochondria and chloroplast
# inspect all high-NegCtr taxa manually
# Remove more potential contaminants
TAX.bac <- TAX.bac[-grep(
  'Eukaryota', TAX.bac$Kingdom),]
TAX.bac <- TAX.bac[-grep(
  'Chloroplast', TAX.bac$Order),] 
TAX.bac <- TAX.bac[-grep(
  'Corynebacteriaceae|Bacillaceae|
  Xanthobacteraceae|Burkholderiaceae|
  Streptococcaceae|Propionibacteriaceae|
  Weeksellaceae|Enterococcaceae', TAX.bac$Family),] 

# match TAX after contaminant removal
ASV.bac <- ASV.bac[row.names(TAX.bac),]

# remove NegCtr columns
ASV.bac <- ASV.bac[, !grepl(
  'neg|NK', names(ASV.bac))]

####################################

## FORMAT TAXONOMY 

# Rename BAC-NAs with last known taxrank + "uc"
k <- ncol(TAX.bac)-1
for (i in 2:k) {
  if (sum(is.na(TAX.bac[, i])) >1) {
    temp <- TAX.bac[is.na(TAX.bac[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX.bac[is.na(TAX.bac[, i]), ] <- temp}
  if (sum(is.na(TAX.bac[, i]))==1) {
    temp <- TAX.bac[is.na(TAX.bac[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX.bac[is.na(TAX.bac[, i]),] <- temp
  }
}
TAX.bac[is.na(TAX.bac[, (k+1)]), (k+1)] <- paste(
  TAX.bac[is.na(TAX.bac[, (k+1)]), k], " uc", sep="")

# shorten/modify names
TAX.bac <- TAX.bac %>%
  mutate(across(everything(),~gsub("Clade_","SAR11 Clade ", .))) %>%
  mutate(across(everything(),~gsub("_clade","", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter_NAC11-7_lineage","NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406)","SAR406", ., fixed=T)))


####################################
## EUKARYOTES / 18S
####################################

ASV.euk <- read.table(
  "./euk_output/eukSeqtab.txt",
  h = T, sep = "\t",
  check.names=F)

TAX.euk <- read.table(
  "./euk_output/eukTax.txt",
  h = T, 
  sep = "\t")

# Rename PR2 taxranks: "Division" to "Phylum"
# OK since taxnames are consistent with Silva
# e.g. Haptophyta: Silva-Phylum, PR2-Division
# enables cross-compatibility with BAC
colnames(TAX.euk)<- c(
  "Kingdom","Supergroup","Phylum","Class",
  "Order","Family","Genus","Species")

# remove Supergroup
TAX.euk$Supergroup <- NULL

#######################################

### CONTAMINANT CHECK

# Define samples PCR-amplified with 25/30/35 cycles
NK25 <- meanDate %>% filter(
  cycles=="25" & locus_tag=="18S") %>% pull(clip_id)
NK30 <- meanDate %>% filter(
  cycles=="30" & locus_tag=="18S") %>% pull(clip_id)
NK35 <- meanDate %>% filter(
  cycles=="35" & locus_tag=="18S") %>% pull(clip_id)

# Define corresponding negative controls
ASV.euk$NK25 <- rowMeans(ASV.euk[,c(
  "NK-28-80_clip","NK-82-126_clip",
  "neg-control-EUK-25c_clip",
  "PS126_NK_25c_EUK_clip")])
ASV.euk$NK30 <- rowMeans(ASV.euk[,c(
  "neg-control-EUK-30c_clip",
  "PS126_NK_30c_EUK_clip")])
ASV.euk$NK35 <- rowMeans(ASV.euk[,c(
  "neg-control-EUK-35c_clip",
  "PS126_NK_35c_EUK_clip")])

# round
ASV.euk$NK25 <- round(ASV.euk$NK25, 0)
ASV.euk$NK30 <- round(ASV.euk$NK30, 0)
ASV.euk$NK35 <- round(ASV.euk$NK35, 0)

# Subtract negative counts 
asv1 = ASV.euk[NK25] - ASV.euk$NK25
asv2 = ASV.euk[NK30] - ASV.euk$NK30
asv3 = ASV.euk[NK35] - ASV.euk$NK35

# Join everything again
ASV.euk <- cbind(asv1, asv2, asv3)

# Set neg to zero
ASV.euk[ASV.euk < 0] <- 0

# Export Animalia/Metazoa to new DF
TAX.meta <- TAX.euk[grep('Metazoa', TAX.euk$Phylum),]
ASV.meta <- ASV.euk[row.names(TAX.meta),]

# Remove Animalia/Metazoa from org table
TAX.euk <- TAX.euk[-grep(
  'Metazoa', TAX.euk$Phylum),]

# match TAX after contaminant removal
ASV.euk <- ASV.euk[row.names(TAX.euk),]

# Remove NegCtr columns 
ASV.euk <- ASV.euk[, !grepl(
  'NK|neg', names(ASV.euk))]
ASV.meta <- ASV.meta[, !grepl(
  'NK|neg', names(ASV.meta))]

####################################

### FORMAT TAXONOMY ###

# Rename EUK-NAs with last known taxrank + "uc"
k <- ncol(TAX.euk)-1
for (i in 2:k) {
  if (sum(is.na(TAX.euk[, i])) >1) {
    temp <- TAX.euk[is.na(TAX.euk[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX.euk[is.na(TAX.euk[, i]), ] <- temp}
  if (sum(is.na(TAX.euk[, i]))==1) {
    temp <- TAX.euk[is.na(TAX.euk[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX.euk[is.na(TAX.euk[, i]),] <- temp
  }
}
TAX.euk[is.na(TAX.euk[, (k+1)]), (k+1)] <- paste(
  TAX.euk[is.na(TAX.euk[, (k+1)]), k], " uc", sep="")

## shorten/modify EUK names 
TAX.euk <- TAX.euk %>%
  mutate(across(everything(),~gsub("Dino-Group-I-Clade-","Dino-I-", .))) %>%
  mutate(across(everything(),~gsub("Dino-Group-II-Clade-","Dino-II-", .))) %>%
  mutate(across(everything(),~gsub("Polar-centric-","", .))) %>%
  mutate(across(everything(),~gsub("Radial-centric-basal-","", .))) %>%
  mutate(across(everything(),~gsub("Chrysophyceae_Clade-","Chrysophyceae ", .))) %>%
  mutate(across(everything(),~gsub("Stephanoecidae_Group_","Stephanoecidae ", .))) %>%
  mutate(across(everything(),~gsub("Pirsonia_Clade_","Pirsonia", .))) %>%
  mutate(across(everything(),~gsub("_X|_XX|_XXX|_XXXX"," uc", .))) 


############################################################################################
###  LOAD + FORMAT ENV-DATA  ###
############################################################################################

## CTD data ##
CTD <- read.table(
  "/AWI_MPI/FRAM/RAS/ampliconTimeseries/metadata/CTD.txt", 
  h=T, sep="\t", stringsAsFactors=F, skipNul=T) #%>%
  #mutate_at(vars(date), as.Date, format = "%Y-%m-%d") 

## NUTRIENTS ##
Nutri <- read.table(
  "/AWI_MPI/FRAM/RAS/ampliconTimeseries/metadata/Nutrients.txt",
  h = T, sep = "\t", check.names=F) 

## STRATIFICATION ##
Strat <- read.table(
  "/AWI_MPI/FRAM/RAS/ampliconTimeseries/metadata/Strat.txt",
  h = T, sep = "\t", check.names=F) %>%
  group_by(date, mooring) %>%
  summarize_at(c("MLD"), mean, na.rm=T) %>%
  ungroup()

## CHLOROPHYLL SATELLITE
chl_sat <- read.table(
  "/AWI_MPI/FRAM/RAS/ampliconTimeseries/metadata/Chl_sat.txt", 
  h=T, sep ="\t",  check.names=F) %>%
  #dplyr::rename(chl_sat=EGC) %>%
  #mutate(mooring="EGC") %>%
  # mutate(date=as.Date(
  #   date, format = "%Y-%m-%d")) %>%
  reshape2::melt(id.vars=c("date")) %>%
  dplyr::rename(
    mooring=variable, chl_sat=value)

####################################

## ICE COVER ##
IceConc <- read.table(
  "/AWI_MPI/FRAM/RAS/ampliconTimeseries/metadata/IceConc.txt",
  h=T, sep = "\t",  check.names=F) %>%
  reshape2::melt() %>%
  dplyr::rename(
    mooring = variable, 
    iceConc = value) %>% 
  left_join(dplyr::select(
    ENV, RAS_id, mooring, date, variable)) %>%
  distinct(date, mooring, .keep_all=T)

# Calculate past ice conditions
# Remove duplicates (16/18S per date)
# fill up NAs with following RAS_id
IceConcPast <- IceConc[!duplicated(
  IceConc$RAS_id, incomparables=NA),] %>%
  fill(RAS_id, .direction='up')

# Calculate past ice
IceConcPast <- IceConcPast %>%
  group_by(RAS_id) %>%
  dplyr::summarize(icePast=mean(iceConc)) %>%
  drop_na() 

# Ice edge distance
IceDist <- read.table(
  "/AWI_MPI/FRAM/RAS/ampliconTimeseries/metadata/IceDist.txt",
  h=T, sep="\t",  check.names=F) %>% 
  reshape2::melt() %>%
  dplyr::rename(
    mooring = variable, 
    iceDist = value) %>%
  mutate_if(is.numeric, round, 2) %>% 
  left_join(dplyr::select(
    ENV, RAS_id, mooring, date, variable)) %>%
  distinct(date, mooring, .keep_all=T)

# Calculate past ice conditions
# remove duplicates (16/18S per date)
# fill up NAs with following RAS_id
IceDistPast <- IceDist[!duplicated(
  IceDist$RAS_id, incomparables=NA),] %>%
  fill(RAS_id, .direction='up')

# Calculate past ice distance
IceDistPast <- IceDistPast %>%
  group_by(RAS_id) %>%
  dplyr::summarize(iceDistPast=mean(iceDist)) %>%
  drop_na() 

# calculate mean ice cover/dist per RAS_id
IceConc <- IceConc %>%
  group_by(RAS_id) %>%
  summarise(iceConc=mean(iceConc)) %>%
  ungroup
IceDist <- IceDist %>%
  group_by(RAS_id) %>%
  summarise(iceDist=mean(iceDist)) %>%
  ungroup
  
##################################

# Join all
# Round numerics
# Convert NaN to NA
ENV <- ENV %>%
  left_join(CTD) %>%
  left_join(Nutri) %>%
  left_join(IceConc) %>%
  left_join(IceDist) %>%
  #left_join(IceConc, by=c("date","mooring_full")) %>% 
  left_join(IceConcPast) %>%
  left_join(IceDistPast) %>%
  left_join(Strat) %>% 
  distinct() %>%
  group_by(RAS_id, locus_tag) %>%
  summarize_if(is.numeric, mean, na.rm=T) %>%
  ungroup() %>% 
  #select_if(~ !all(is.na(.))) %>%
  left_join(meanDate) %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .))

# Correct timecourse
ENV$month_full <- factor(ENV$month_full, levels=c(
  "Jul-16","Aug-16","Sep-16","Oct-16",
  "Nov-16","Dec-16","Jan-17","Feb-17",
  "Mar-17","Apr-17","May-17","Jun-17",
  "Jul-17","Aug-17","Sep-17","Oct-17",
  "Nov-17","Dec-17","Jan-18","Feb-18",
  "Mar-18","Apr-18","May-18","Jun-18",
  "Jul-18","Aug-18","Sep-18","Oct-18",
  "Nov-18","Dec-18","Jan-19","Feb-19",
  "Mar-19","Apr-19","May-19","Jun-19",
  "Jul-19","Aug-19","Sep-19","Oct-19",
  "Nov-19","Dec-19","Jan-20","Feb-20",
  "Mar-20","Apr-20","May-20","Jun-20",
  "Jul-20","Aug-20","Sep-20"))

# separate EUK/BAC
ENV.bac <- as.data.frame(ENV) %>%
  filter(locus_tag=="16S")
ENV.euk <- as.data.frame(ENV) %>%
  filter(locus_tag=="18S")


############################################################################################
  ### SORT + ALIGN WITH ASVs ###
############################################################################################

# Sort data 
ASV.bac <- ASV.bac[,mixedsort(names(ASV.bac))]
ENV.bac <- ENV.bac[mixedorder(ENV.bac$clip_id),]

# Confirm consistency of data  
ASV.bac <- ASV.bac[,c((match(
  ENV.bac$clip_id, colnames(ASV.bac))))]

# Rename 
colnames(ASV.bac) = ENV.bac$RAS_id

##########################

# Sort data 
ASV.euk <- ASV.euk[,mixedsort(names(ASV.euk))]
ENV.euk <- ENV.euk[mixedorder(ENV.euk$clip_id),]

# Confirm consistency of data  
ASV.euk <- ASV.euk[,c((match(
  ENV.euk$clip_id, colnames(ASV.euk))))]

# Rename 
colnames(ASV.euk) = ENV.euk$RAS_id


############################################################################################
### EXPORT MASTER-TABLES ###
############################################################################################

## 16S
write.table(
  ASV.bac, file="bacASV.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
write.table(
  TAX.bac, file="bacTAX.txt",
  sep="\t", row.names=T, col.names=T, quote=F)

## 18S
write.table(
  ASV.euk, file="eukASV.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
write.table(
  TAX.euk, file="eukTAX.txt",
  sep="\t", row.names=T, col.names=T, quote=F)

## 18S -- metazoans
write.table(
  ASV.meta, file="metASV.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
write.table(
  TAX.meta, file="metTAX.txt",
  sep="\t", row.names=F, col.names=T, quote=F)

# Metadata
write.table(
  ENV, file="_metadata.txt",
  sep="\t", row.names=F, col.names=T, quote=F)
write.table(
  ENV.bac, file="bacMeta.txt",
  sep="\t", row.names=F, col.names=T, quote=F)
write.table(
  ENV.euk, file="eukMeta.txt",
  sep="\t", row.names=F, col.names=T, quote=F)

############################################################################################

# remove temp-data
rm(asv1, asv2, asv3, NK25, NK30, NK35,
   temp, i, j, k, subset, asvB, asvE)

# save Rdata
save.image("_RAS.Rdata")
