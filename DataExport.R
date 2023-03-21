#########################################################
### FRAM-RAS -- AMPLICON TIME-SERIES  ###
#########################################################

### EXPORT FOR SPECIFIC STUDIES 
#Script describes how to subset the data for specific studies 

library(dplyr)
library(tidyr)


#########################################################
## EGC XX: doi() ##
#########################################################

ASV.bac %>%
  dplyr::select(c(contains("EGC"))) %>%
  write.table(
    file="../MetaG-Ice/Rstats/FRAM_RAS_EGC_ASV_raw.txt",
    sep="\t", row.names=T, col.names=T, quote=F)
TAX.bac %>% write.table(
  file="../MetaG-Ice/Rstats/FRAM_RAS_EGC_ASV_taxa.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
ENV %>%
  filter(mooring=="EGC" & locus_tag=="16S") %>%
  write.table(
    file="../MetaG-Ice/Rstats/FRAM_RAS_EGC_ASV_meta.txt",
    sep="\t", row.names=F, col.names=T, quote=F)

######################################################
######################################################

## WSC / F4 XX: doi() ##

ENV %>%
  filter(mooring=="F4") %>% write.table(
    file="metadata_WSC.txt",
    sep="\t", row.names=F, col.names=T, quote=F)

#######################

ASV.bac %>%
  dplyr::select(c(contains("F4"))) %>% write.table(
    file="../MetaG-WSC/Rstats/bacASV.txt",
    sep="\t", row.names=T, col.names=T, quote=F)
TAX.bac %>% write.table(
  file="../MetaG-WSC/Rstats/bacTAX.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
ENV.bac %>%
  filter(mooring=="F4") %>%
  write.table(
    file="../MetaG-WSC/Rstats/bacMeta.txt",
    sep="\t", row.names=F, col.names=T, quote=F)

#######################

ASV.euk %>%
  dplyr::select(c(contains("F4"))) %>%
  write.table(
    file="../MetaG-WSC/Rstats/eukASV.txt",
    sep="\t", row.names=T, col.names=T, quote=F)
TAX.euk %>% write.table(
  file="../MetaG-WSC/Rstats/eukTAX.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
ENV.euk %>%
  filter(mooring=="F4") %>%
  write.table(
    file="../MetaG-WSC/Rstats/eukMeta.txt",
    sep="\t", row.names=F, col.names=T, quote=F)


#########################################################
## HG-IV+F4 2017-18: doi() ##
#########################################################

subset <- ENV %>%
  filter(locus_tag=="16S" & mooring_full %in% c(
    "F4-S-2","HG-IV-S-2")) %>%
  pull(RAS_id)

ASV.bac %>% dplyr::select(subset) %>% write.table(
  file="../Ellen_1718/bacASV.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
TAX.bac %>% write.table(
  file="../Ellen_1718/bacTAX.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
ENV.bac %>%
  filter(mooring_full %in% c("F4-S-2","HG-IV-S-2")) %>% write.table(
    file="../Ellen_1718/bacMeta.txt",
    sep="\t", row.names=F, col.names=T, quote=F)

#######################

subset <- ENV %>%
  filter(locus_tag=="18S" & mooring_full %in% c(
    "F4-S-2","HG-IV-S-2")) %>%
  pull(RAS_id)

ASV.euk %>% dplyr::select(subset) %>% write.table(
  file="../Ellen_1718/eukASV.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
TAX.euk %>% write.table(
  file="../Ellen_1718/eukTAX.txt",
  sep="\t", row.names=T, col.names=T, quote=F)
ENV.euk %>%
  filter(mooring_full %in% c("F4-S-2","HG-IV-S-2")) %>% write.table(
    file="../Ellen_1718/eukMeta.txt",
    sep="\t", row.names=F, col.names=T, quote=F)


#########################################################
## polarDNAexplorer ##
#########################################################

# Remove unnecessary columns
# Rename variables for clarity
# Sort correctly
ENV %>% dplyr::select(-c(
  "date1","date2","date3","date4","sig",
  "CO2","pH","clip_id","cycles","pcr_primer",
  "target_subfragment","DNA_ng_uL",
  "miseq_id","PangaeaEvent_Deployment",
  "PangaeaEvent_Recovery","read_fwd","read_rev",
  "read_fwd_run2","read_rev_run2")) %>% 
  arrange(locus_tag, mooring_full, date) %>%
  dplyr::rename(
    `rRNA gene` = locus_tag,
    `Sea ice cover (%)` = iceConc,
    `Sea ice cover (%; past)` = icePast,
    `Distance to ice edge (km)` = iceDist,
    `Distance to ice edge (km; past)` = iceDistPast,
    `Proportion of Atlantic Water (%)` = AW_frac,
    `Proportion of Polar Water (%)` = PW_frac,
    `Mixed layer depth (m)` = MLD,
    `Water temperature (°C)` = temp,
    `Salinity (psu)` = sal,
    `Daylight (h)` = daylight,
    `Degree Latitude` = lat,
    `Degree Longitude` = lon,
    `Oxygen concentration (mg/L)` = O2_conc,
    `Oxygen saturation (%)` = O2_sat,
    `Sensor chlorophyll-a` = chl_sens,
    `Nitrate + nitrite (µM)` = NO3_NO2,
    `Phosphate (µM)` = PO4,
    `Nitrate (µM)` = NO2,
    `Silicate (µM)` = SiO4,
    `Sampling depth (m)` = depth,
    `Sampling date` = date,
    `Mooring ID` = mooring_full,
    Year = year) %>% write.table(
      file="../../polarDNAexplorer/metadata.txt",
      sep="\t", row.names=F, col.names=T, quote=F)

##################################

# remove low-abundance ASVs (as done in all specific studies)
# Subset TAX accordingly

asvB <- ASV.bac %>% 
  filter(rowSums(.>= 3) >= 3) 
write.table(
  asvB, file="../../polarDNAexplorer/bacASV.txt",
  sep="\t", row.names=T, col.names=T, quote=F)

TAX.bac %>% 
  filter(rownames(.) %in% rownames(asvB)) %>% 
  write.table(
    file="../../polarDNAexplorer/bacTAX.txt",
    sep="\t", row.names=T, col.names=T, quote=F)

asvE <- ASV.euk %>% 
  filter(rowSums(.>= 3) >= 3) 
write.table(
  asvE, file="../../polarDNAexplorer/eukASV.txt",
  sep="\t", row.names=T, col.names=T, quote=F)

TAX.euk %>% filter(rownames(.) %in% rownames(asvE)) %>% 
  write.table(
    file="../../polarDNAexplorer/eukTAX.txt",
    sep="\t", row.names=T, col.names=T, quote=F)


#########################################################
# SPLIT VS SEPARATE ASVs (doi)
#########################################################

ENV.bac %>%
  filter(mooring_full %in% c("F4-S-1","EGC-5")) %>% write.table(
    file="../../../collaborations/SplitPrimer/bacMeta.txt",
    sep="\t", row.names=F, col.names=T, quote=F)
ENV.euk %>%
  filter(mooring_full %in% c("F4-S-1","EGC-5")) %>% write.table(
    file="../../../collaborations/SplitPrimer/eukMeta.txt",
    sep="\t", row.names=F, col.names=T, quote=F)

