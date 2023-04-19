# Global M-Value Correlations

# Load packages
library(dplyr)
library(missMethyl)
library(ENmix)
library(limma)
library(ggplot2)
library(DMRcate)
library(readr)
library(piggyback)
library(tidyverse)
# library(lme4)

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

###### DMP and DMR Identification in All Transgender Youth on GAH ######
# Load objects from "preprocessing_and_normalization.R" code
# G954 metadata, cell type, and methylation age
G954_mage_ageAcc <- readRDS("G954_mage_ageAcc.RDS")
# G954 RGChannelSetExtended
rg <- readRDS("G954_EPIC_Data_0217rgsetext.RDS")
# G954 beta values
rg_betas <- readRDS("G954_EPIC_Data_0217rg_betas.RDS")

# Winsorize beta values (adapted from Benjamin Jacobson's thesis)
rg_betas[rg_betas > .99] <- .99
rg_betas[rg_betas < .01] <- .01

# Convert beta values to M-values using ENmix
Mval<-B2M(rg_betas)

# Remove metadata columns for samples not of interest (i.e., paired data 
# unusable after ENmix quality control)
G954<-G954_mage_ageAcc %>% 
    subset(., !(Sample_Name %in% c("206839580077_R08C01",# 10007_2
                                   "206839580077_R07C01")))# 10007
# Select only metadata rows for people on GAH
G954_GAH<- G954[G954$GAH %in% 1,]

# Create metadata objects for baseline
G954_GAH_baseline <- G954_GAH[G954_GAH$Visit %in% 1,]
# Create metadata objects for Visit 2
G954_GAH_2 <- G954_GAH[G954_GAH$Visit %in% 2,]
# Create metadata object for Visit 3
G954_GAH_3<-G954_GAH[G954_GAH$Visit %in% 3,]

# Create vector of Sample_Names from baseline 
string_use_baseline <- G954_GAH_baseline[["Sample_Name"]]
# Create vector of Sample_Names from Visit 2 
string_use_2 <- G954_GAH_2[["Sample_Name"]]
# Create vector of Sample_Names from Visit 3
string_use_3 <- G954_GAH_3[,"Sample_Name"]

# Select M-value columns for cross-sectional (baseline) samples 
Mval_baseline<-Mval[,(colnames(Mval) %in% string_use_baseline)] 
# Select M-value columns for cross-sectional (Visit 2) samples 
Mval_Visit2<-Mval[,(colnames(Mval) %in% string_use_2)] 
# Select M-value columns for cross-sectional (Visit 3) samples 
Mval_Visit3<-Mval[,(colnames(Mval) %in% string_use_3)] 

# Load manifest (from Ben's thesis)
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/zhou")
man <- readRDS("EPIC.hg19.manifest.rds")
man <- as(man, "data.frame")
man$Probe_ID <- rownames(man)
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

# Exclude bad probes and sex chromosomes (from Ben's thesis)
manSub <- man[!man$MASK_general,]
manSub <- manSub[manSub$seqnames != "chrX",]
manSub <- manSub[manSub$seqnames != "chrY",]
manSub <- manSub[manSub$seqnames != "chrM",]

# Exclude all bad probes, SNPs, and sex chromosomes from M-value data (from Ben's 
# thesis)
# Cross-sectional (baseline)
Mval_baseline_clean<-Mval_baseline[rownames(Mval_baseline) %in% manSub$Probe_ID,]
# Cross-sectional (Visit 2)
Mval_Visit2_clean<-Mval_Visit2[rownames(Mval_Visit2) %in% manSub$Probe_ID,]
# Cross-sectional (Visit 3)
Mval_Visit3_clean<-Mval_Visit3[rownames(Mval_Visit3) %in% manSub$Probe_ID,]


# Remove probes that don't change: cross-sectional (baseline)
std<-apply(Mval_baseline_clean,1,sd)
cutoff<-summary(std)[2]/2
keep<-std>cutoff
Mval_baseline2<-Mval_baseline_clean[keep,]

# Remove probes that don't change: cross-sectional (Visit 2)
std<-apply(Mval_Visit2_clean,1,sd)
cutoff<-summary(std)[2]/2
keep<-std>cutoff
Mval_Visit22<-Mval_Visit2_clean[keep,]

# Remove probes that don't change: cross-sectional (Visit 3)
std<-apply(Mval_Visit3_clean,1,sd)
cutoff<-summary(std)[2]/2
keep<-std>cutoff
Mval_Visit33<-Mval_Visit3_clean[keep,]

# Question 1a
# Does Visit 1 PSS10 stress level correlate with Visit 1 DNAm in those on GAH?
# n = 12
# Test it:
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
# Load data
PSSdata<-read.csv("PSSdata.csv")
# Create dataframe with sample id and baseline PSS values
baseline_PSS <- PSSdata[PSSdata$Visit==1,c("Subject.ID","PSS")]
names(baseline_PSS)[2]<-"baseline_PSS"
# Merge Visit 1 metadata with baseline_PSS
pssmeta1<-merge(baseline_PSS, G954_GAH_baseline, by="Subject.ID") %>% arrange(Sample_Name)
# Create string of Sample_Names from pssmeta2
string_pssmeta1<-pssmeta1[["Sample_Name"]]
# Subset Mval_Visit2 to only include columns for those in pssmeta
Mval_baseline_pss<- Mval_baseline2[,(colnames(Mval_baseline2) %in% string_pssmeta1)] 
# Transpose Mval data
Mval_baseline_psst<-t(Mval_baseline_pss)
# Make dataframe
Mval_baseline_psst<-as.data.frame(Mval_baseline_psst)
# Make new column with row means of methylation change
Mval_baseline_psst$rowmeans <- rowMeans(Mval_baseline_psst,
                                        na.rm=TRUE)
# Make column for rownames
Mval_baseline_psst$Sample_Name<-rownames(Mval_baseline_psst)
# Merge Mval data with PSS data
Mval1pss<-merge(Mval_baseline_psst,pssmeta1, by="Sample_Name")
# Create scatterplot
ggplot(Mval1pss, aes(x = rowmeans, y = baseline_PSS)) +
    geom_point()
# Save csv
Mval1pss.csv<-Mval1pss %>% select(Sample_Name,rowmeans,baseline_PSS)
write.csv(Mval1pss.csv,file="Mval1pss.csv")


# Question 1b
# Does Visit 2 PSS10 stress level correlate with Visit 2 DNAm in those on GAH?
# n = 13
# Test it:
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
# Load data
PSSdata<-read.csv("PSSdata.csv")
# Create dataframe with sample id and Visit 2 PSS values
Visit2_PSS <- PSSdata[PSSdata$Visit==2,c("Subject.ID","PSS")]
names(Visit2_PSS)[2]<-"Visit2_PSS"
# Merge Visit 2 metadata with Visit2_PSS
pssmeta2<-merge(Visit2_PSS, G954_GAH_2, by="Subject.ID") %>% arrange(Sample_Name)
# Create string of Sample_Names from pssmeta2
string_pssmeta2<-pssmeta2[["Sample_Name"]]
# Subset Mval_Visit2 to only include columns for those in pssmeta
Mval_Visit2_pss<- Mval_Visit22[,(colnames(Mval_Visit22) %in% string_pssmeta2)] 
# Transpose Mval data
Mval_Visit2_psst<-t(Mval_Visit2_pss)
# Make dataframe
Mval_Visit2_psst<-as.data.frame(Mval_Visit2_psst)
# Make new column with row means of methylation change
Mval_Visit2_psst$rowmeans <- rowMeans(Mval_Visit2_psst,
                                      na.rm=TRUE)
# Make column for rownames
Mval_Visit2_psst$Sample_Name<-rownames(Mval_Visit2_psst)
# Merge Mval data with PSS data
Mval2pss<-merge(Mval_Visit2_psst,pssmeta2, by="Sample_Name")
# Create scatterplot
ggplot(Mval2pss, aes(x = rowmeans, y = Visit2_PSS)) +
    geom_point()
# Save csv
Mval2pss.csv<-Mval2pss %>% select(Sample_Name,rowmeans,Visit2_PSS)
write.csv(Mval2pss.csv,file="Mval2pss.csv")

# Question 1c
# Does Visit 3 PSS10 stress level correlate with Visit 3 DNAm in those on GAH?
# n = 3
# Test it:
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
# Load data
PSSdata<-read.csv("PSSdata.csv")
# Create dataframe with sample id and Visit 2 PSS values
Visit3_PSS <- PSSdata[PSSdata$Visit==3,c("Subject.ID","PSS")]
names(Visit3_PSS)[2]<-"Visit3_PSS"
# Merge Visit 2 metadata with Visit2_PSS
pssmeta3<-merge(Visit3_PSS, G954_GAH_3, by="Subject.ID") %>% arrange(Sample_Name)
# Create string of Sample_Names from pssmeta2
string_pssmeta3<-pssmeta3[["Sample_Name"]]
# Subset Mval_Visit2 to only include columns for those in pssmeta
Mval_Visit3_pss<- Mval_Visit33[,(colnames(Mval_Visit33) %in% string_pssmeta3)] 
# Transpose Mval data
Mval_Visit3_psst<-t(Mval_Visit3_pss)
# Make dataframe
Mval_Visit3_psst<-as.data.frame(Mval_Visit3_psst)
# Make new column with row means of methylation change
Mval_Visit3_psst$rowmeans <- rowMeans(Mval_Visit3_psst,
                                      na.rm=TRUE)
# Make column for rownames
Mval_Visit3_psst$Sample_Name<-rownames(Mval_Visit3_psst)
# Merge Mval data with PSS data
Mval3pss<-merge(Mval_Visit3_psst,pssmeta3, by="Sample_Name")
# Create scatterplot
ggplot(Mval3pss, aes(x = rowmeans, y = Visit3_PSS)) +
    geom_point()
# Save csv
Mval3pss.csv<-Mval3pss %>% select(Sample_Name,rowmeans,Visit3_PSS)
write.csv(Mval3pss.csv,file="Mval3pss.csv")


# Question 2a
# Does Visit 1 PROMIS anxiety level correlate with Visit 1 DNAm in those on GAH?
# n = 12
# Test it:
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
# Load PROMIS data
PROMISdata<-read.csv("PROMISdata.csv")
# Create dataframe with sample id and baseline PROMIS Anxiety values
baseline_PROMISanx <- PROMISdata[PROMISdata$Visit==1,
                                 c("Subject.ID",
                                   "PROMIS_ped_anx_tscore",
                                   "PROMIS_adult_anx_tscore")]
baseline_PROMISanx %>% unite(col="baseline_PROMISanx",PROMIS_ped_anx_tscore,
                             PROMIS_adult_anx_tscore,na.rm=T) ->
    baseline_PROMISanx
# Make column numeric
baseline_PROMISanx$baseline_PROMISanx=as.numeric(baseline_PROMISanx$baseline_PROMISanx)
# Merge Visit 1 metadata with PROMISanx_1_2
promismeta1<-merge(baseline_PROMISanx, G954_GAH_baseline, by="Subject.ID") %>% 
    arrange(Sample_Name)
# Create string of Sample_Names from promismeta
string_promismeta1<-promismeta1[["Sample_Name"]]
# Subset Mval_Visit2 to only include columns for those in promismeta2
Mval_baseline_promis<- Mval_baseline2[,(colnames(Mval_baseline2) %in% string_promismeta1)] 
# Transpose Mval data
Mval_baseline_promist<-t(Mval_baseline_promis)
# Make dataframe
Mval_baseline_promist<-as.data.frame(Mval_baseline_promist)
# Make new column with row means of methylation change
Mval_baseline_promist$rowmeans <- rowMeans(Mval_baseline_promist,
                                           na.rm=TRUE)
# Make column for rownames
Mval_baseline_promist$Sample_Name<-rownames(Mval_baseline_promist)
# Merge Mval data with PROMIS data
Mval1promis<-merge(Mval_baseline_promist,promismeta1, by="Sample_Name")
# Create scatterplot
ggplot(Mval1promis, aes(x = rowmeans, y = baseline_PROMISanx)) +
    geom_point()
# Save csv
Mval1promis.csv<-Mval1promis %>% select(Sample_Name,rowmeans,baseline_PROMISanx)
write.csv(Mval1promis.csv,file="Mval1promis.csv")

# Question 2b
# Does Visit 2 PROMIS anxiety level correlate with Visit 2 DNAm in those on GAH?
# n = 13
# Test it:
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
# Load PROMIS data
PROMISdata<-read.csv("PROMISdata.csv")
# Create dataframe with sample id and Visit 2 PROMIS Anxiety values
Visit2_PROMISanx <- PROMISdata[PROMISdata$Visit==2,
                               c("Subject.ID",
                                 "PROMIS_ped_anx_tscore",
                                 "PROMIS_adult_anx_tscore")]
Visit2_PROMISanx %>% unite(col="Visit2_PROMISanx",PROMIS_ped_anx_tscore,
                           PROMIS_adult_anx_tscore,na.rm=T) ->Visit2_PROMISanx
# Make column numeric
Visit2_PROMISanx$Visit2_PROMISanx=as.numeric(Visit2_PROMISanx$Visit2_PROMISanx)
# Merge Visit 2 metadata with PROMISanx_1_2
promismeta2<-merge(Visit2_PROMISanx, G954_GAH_2, by="Subject.ID") %>% 
    arrange(Sample_Name)
# Create string of Sample_Names from promismeta
string_promismeta2<-promismeta2[["Sample_Name"]]
# Subset Mval_Visit2 to only include columns for those in promismeta2
Mval_Visit2_promis<- Mval_Visit22[,(colnames(Mval_Visit22) %in% string_promismeta2)] 
# Transpose Mval data
Mval_Visit2_promist<-t(Mval_Visit2_promis)
# Make dataframe
Mval_Visit2_promist<-as.data.frame(Mval_Visit2_promist)
# Make new column with row means of methylation change
Mval_Visit2_promist$rowmeans <- rowMeans(Mval_Visit2_promist,
                                         na.rm=TRUE)
# Make column for rownames
Mval_Visit2_promist$Sample_Name<-rownames(Mval_Visit2_promist)
# Merge Mval data with PROMIS data
Mval2promis<-merge(Mval_Visit2_promist,promismeta2, by="Sample_Name")
# Create scatterplot
ggplot(Mval2promis, aes(x = rowmeans, y = Visit2_PROMISanx)) +
    geom_point()
# Save csv
Mval2promis.csv<-Mval2promis %>% select(Sample_Name,rowmeans,Visit2_PROMISanx)
write.csv(Mval2promis.csv,file="Mval2promis.csv")

# Question 2c
# Does Visit 3 PROMIS anxiety level correlate with Visit 3 DNAm in those on GAH?
# n = 3
# Test it:
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
# Load PROMIS data
PROMISdata<-read.csv("PROMISdata.csv")
# Create dataframe with sample id and Visit 2 PROMIS Anxiety values
Visit3_PROMISanx <- PROMISdata[PROMISdata$Visit==3,
                               c("Subject.ID",
                                 "PROMIS_ped_anx_tscore",
                                 "PROMIS_adult_anx_tscore")]
Visit3_PROMISanx %>% unite(col="Visit2_PROMISanx",PROMIS_ped_anx_tscore,
                           PROMIS_adult_anx_tscore,na.rm=T) ->Visit3_PROMISanx
# Merge Visit 2 metadata with PROMISanx_1_2
promismeta3<-merge(Visit3_PROMISanx, G954_GAH_3, by="Subject.ID") %>% 
    arrange(Sample_Name)
# Create string of Sample_Names from promismeta
string_promismeta3<-promismeta3[["Sample_Name"]]
# Subset Mval_Visit2 to only include columns for those in promismeta2
Mval_Visit3_promis<- Mval_Visit33[,(colnames(Mval_Visit33) %in% string_promismeta3)] 
# Transpose Mval data
Mval_Visit3_promist<-t(Mval_Visit3_promis)
# Make dataframe
Mval_Visit3_promist<-as.data.frame(Mval_Visit3_promist)
# Make new column with row means of methylation change
Mval_Visit3_promist$rowmeans <- rowMeans(Mval_Visit3_promist,
                                         na.rm=TRUE)
# Make column for rownames
Mval_Visit3_promist$Sample_Name<-rownames(Mval_Visit3_promist)