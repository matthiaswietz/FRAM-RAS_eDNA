
############################################################################################
   ###  RAS -- AMPLICON ANALYSIS  ###
############################################################################################

# This script: rarefaction and alpha-diversity 

library(iNEXT)
library(olsrr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tibble)
library(gtools)
#load("_RAS_DataLoad.Rdata")
load("iNEXT_bac.RData")
load("iNEXT_euk.RData")

############################################################################################
###  RAREFACTION AND COVERAGE  ###
############################################################################################

#iNEXT.bac <- otu_table(
# ASV.bac, taxa_are_rows=F)
iNEXT.bac <- iNEXT(
  ASV.bac, datatype="abundance", 
  q=c(0), conf=0.95, nboot=100)

iNEXT.euk <- iNEXT(
  ASV.euk, datatype="abundance", 
  q=c(0), conf=0.95, nboot=100)


###################################

## RAREFACTION ##

rarefac.bac <- fortify(iNEXT.bac, type=1) %>%
  right_join(ENV.bac, by=c("site"="clip_id"))

rarefac.line.bac <- rarefac.bac[which(
  rarefac.bac$method != "observed"),]
rarefac.line.bac$method <- factor(
  rarefac.line.bac$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

rarefac.euk <- fortify(iNEXT.euk, type=1) %>%
  right_join(ENV.euk, by=c("site"="clip_id"))

rarefac.line.euk <- rarefac.euk[which(
  rarefac.euk$method != "observed"),]
rarefac.line.euk$method <- factor(
  rarefac.line.euk$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

###################################

## COVERAGE ##

cover.bac <- fortify(iNEXT.bac, type=2) %>%
  right_join(ENV.bac, by=c("site"="clip_id"))

cover.line.bac <- cover.bac [which(
  cover.bac$method != "observed"),]
cover.line.bac$method <- factor(
  cover.line.bac$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

cover.euk <- fortify(iNEXT.euk, type=2) %>%
  right_join(ENV.euk, by=c("site"="clip_id"))

cover.line.euk <- cover.euk [which(
  cover.euk$method != "observed"),]
cover.line.euk$method <- factor(
  cover.line.euk$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

###################################

## COMBINE + PLOT

coverage <- rbind(cover.bac,cover.euk)
cover.line <- rbind(cover.line.bac,cover.line.euk)
rarefaction <- rbind(rarefac.bac,rarefac.euk)
rarefac.line <- rbind(rarefac.line.bac,rarefac.line.euk)

ggplot(coverage, aes(
  x=x, y=y, colour=site))+ 
geom_line(
  aes(linetype=method), 
  lwd = 0.5, data=cover.line) +
scale_colour_discrete(guide ="none") +
scale_x_continuous(
  limits = c(0,1e+5)) +
scale_y_continuous(
  breaks = seq(0.9,1,0.05), 
  limits = c(0.9,1)) +
facet_grid(locus_tag~mooring_full) +
labs(x="Sample size", y="Sample coverage") +
theme_bw(base_size = 12) + 
theme(legend.position="bottom")

ggplot(rarefaction, aes(
  x=x, y=y, colour=site)) +
geom_line(aes(
  linetype=method), 
  lwd=0.5, data=rarefac.line) +
scale_colour_discrete(guide = "none") +
scale_x_continuous(limits = c(0,1e+5)) +
facet_grid(locus_tag~mooring_full) +
labs(x="Sample size", y="Species richness") +
theme_bw() + 
theme(legend.position="bottom")


###########################################################################################
   ###  REFORMATTING ###
############################################################################################

richness <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Species richness",] %>%
  arrange(Site) %>%
  right_join(ENV.bac, by=c("Site"="clip_id"))

simpson <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Site) %>%
  right_join(ENV.bac, by=c("Site"="clip_id"))

div.bac <- data.frame(
  RAS_id = ENV.bac$RAS_id,
  Site = ENV.bac$clip_id, 
  locus_tag = ENV.bac$locus_tag,
  richness = richness$Observed,
  simpson = simpson$Observed) 

###################################

richness <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Species richness",] %>%
  arrange(Site) %>%
  right_join(ENV.euk, by=c("Site"="clip_id"))

simpson <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Site) %>%
  right_join(ENV.euk, by=c("Site"="clip_id"))

div.euk <- data.frame(
  RAS_id = ENV.euk$RAS_id,
  Site = ENV.euk$clip_id, 
  locus_tag = ENV.euk$locus_tag,
  richness = richness$Observed,
  simpson = simpson$Observed) 

###################################

###  MERGE ###

AlphaDiv <- rbind(
  div.bac, div.euk) %>%
  left_join(ENV)


###########################################################################################
###   export for individual studies ###
###########################################################################################

# F4 2016-20 (doi)

AlphaDiv %>%
  filter(mooring=="F4") %>%
  write.table(
    file="../MetaG-WSC/Rstats/AlphaDiv.txt",
    sep="\t", row.names=F, col.names=T, quote=F)
  

###############################################################

# remove temp-data
rm(richness, simpson, 
   coverage, cover.bac, cover.euk, 
   rarefaction, rarefac.bac, rarefac.euk,
   div.bac, div.euk)

rm(list = ls(pattern =
 ".point.*|.line.*"))

# Save data
save(iNEXT.bac, file="iNEXT_bac.RData")
save(iNEXT.euk, file="iNEXT_euk.RData")


########

#load("/AWI_MPI/FRAM/RAS/ampliconTimeseries/iNEXT_bac.Rdata")
#load("/AWI_MPI/FRAM/RAS/ampliconTimeseries/iNEXT_euk.Rdata")
