# Age Acceleration Calculations Using Calculations from ENmix and Lisa Schneper

# Load packages
library(dplyr)

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

## Load epigenetic age calculations
# GrimAge data calculated by Lisa Schneper
G954_grimage_outputLS<-readRDS("G954_grimage_outputLS.RDS")
# Horvath, Hannum, and PhenoAge data calculated in ENmix
G954_mage_ageAcc<-readRDS("G954_mage_ageAcc.RDS")

# Attach metadata for later use
baseDir <- file.path("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-20-23")
G954sheet <- read.metharray.sheet(base=baseDir,pattern = "G954_Sample_Sheet.csv")
rownames(G954sheet)=G954sheet$Sample_Name
colData(rg)=as(G954sheet[colnames(rg),],"DataFrame")
G954sheet %>% as.data.frame() %>% 
    arrange(Sample_Name)->G954sheet

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
saveRDS(G954sheet, file = "G954sheet.RDS")

# Calculate Horvath epigenetic age acceleration with ENmix data
EEAAHorvath.ENmix=residuals(lm(mAge_Hovath ~ Age +
                                   Bcell +
                                   CD4T +
                                   CD8T +
                                   Mono +
                                   Neu +
                                   NK, 
                               data=G954_mage_ageAcc, 
                               na.action=na.exclude ))
# Create dataframe
EEAAHorvath.ENmix<-as.data.frame(EEAAHorvath.ENmix)
# Add sample ids
EEAAHorvath.ENmix$Sample_Name<-G954_mage_ageAcc$Sample_Name


# Calculate Hannum epigenetic age acceleration with ENmix data
EEAAHannum.ENmix=residuals(lm(mAge_Hannum ~ Age +
                                  Bcell +
                                  CD4T +
                                  CD8T +
                                  Mono +
                                  Neu +
                                  NK, 
                              data=G954_mage_ageAcc, 
                              na.action=na.exclude ))
# Create dataframe
EEAAHannum.ENmix<-as.data.frame(EEAAHannum.ENmix)
# Add sample ids
EEAAHannum.ENmix$Sample_Name<-G954_mage_ageAcc$Sample_Name


# Calculate PhenoAge epigenetic age acceleration with ENmix data
EEAAPhenoAge.ENmix=residuals(lm(PhenoAge ~ Age +
                                    Bcell +
                                    CD4T +
                                    CD8T +
                                    Mono +
                                    Neu +
                                    NK, 
                                data=G954_mage_ageAcc, 
                                na.action=na.exclude ))
# Create dataframe
EEAAPhenoAge.ENmix<-as.data.frame(EEAAPhenoAge.ENmix)
# Add sample ids
EEAAPhenoAge.ENmix$Sample_Name<-G954_mage_ageAcc$Sample_Name


# Calculate GrimAge epigenetic age acceleration with data from Lisa Schneper
EEAAGrimAge= residuals(lm(DNAmGrimAge~Age +
                              PlasmaBlast +
                              CD8pCD28nCD45RAn +
                              CD8_naive, 
                          data=G954_grimage_outputLS,
                          na.action=na.exclude))
# Create dataframe
EEAAGrimAge<-as.data.frame(EEAAGrimAge)
# Add sample ids
EEAAGrimAge$SampleID<-G954_mage_ageAcc$SampleID
# Rename column for downstream analysis
EEAAGrimAge$Sample_Name<-EEAAGrimAge$SampleID

# Merge all four epigenetic age accelerations
EEAA.ENmix<-merge(EEAAHorvath.ENmix,EEAAHannum.ENmix,by.x="Sample_Name",
                  by.y="Sample_Name")
EEAA.ENmix<-merge(EEAA.ENmix,EEAAPhenoAge.ENmix,by.x="Sample_Name",
                  by.y="Sample_Name")
EEAA.ENmix<-merge(EEAA.ENmix,EEAAGrimAge,by.x="Sample_Name",
                  by.y="Sample_Name")
# Remove SampleID column
EEAA.ENmix<-EEAA.ENmix[,colnames(EEAA.ENmix)!="SampleID"]
# Merge epigenetic age acceleration data with relevant metadata columns
Subject.id.Visit<-G954sheet %>% select(Sample_Name,Subject.ID,Visit)
EEAA.ENmix<-merge(EEAA.ENmix,Subject.id.Visit,by.x="Sample_Name",
                  by.y="Sample_Name")

# Move metadata columns to front
EEAA.ENmix<-EEAA.ENmix %>% relocate(Sample_Name,Subject.ID, Visit)
View(EEAA.ENmix)
# Save as csv
write.csv(EEAA.ENmix,file="EEAA.ENmix.csv")


####### Check Similarity Between Hannum Ages Calculated by Lisa vs ENmix #####
# Calculate Hannum epigenetic age acceleration with Lisa's data
EEAAHannumAge= residuals(lm(DNAmHannumAge~Age 
                            +PlasmaBlast 
                            +CD8pCD28nCD45RAn 
                            +CD8_naive, 
                            data=G954_grimage_outputLS,
                            na.action=na.exclude))
# Create dataframe
EEAAHannumAge<-as.data.frame(EEAAHannumAge)
# Add sample ids
EEAAHannumAge$SampleID<-G954_mage_ageAcc$SampleID
# Rename sample column
EEAAHannumAge$Sample_Name<-EEAAHannumAge$SampleID

# Create dataframe with sample id and baseline epigenetic age acceleration ENmix 
# values
baseline_HannumENmix <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Sample_Name",
                                                         "EEAAHannum.ENmix")]
names(baseline_HannumENmix)[2]<-"baseline_EEAAHannum.ENmix"
# Create dataframe with sample id and baseline epigenetic age acceleration Lisa
# values
baseline_HannumLisa <- EEAAHannumAge[,c("Sample_Name","EEAAHannumAge")]
names(baseline_HannumLisa)[2]<-"baseline_EEAAHannumAge.Lisa"
# Merge into one dataframe with 3 columns
EEAA_Hannum_1_1 <- merge(baseline_HannumENmix, baseline_HannumLisa, 
                         by.x="Sample_Name",by.y="Sample_Name")
# Delete rows by name (delete those only on blockers: 10002, 10017, 10018)
EEAA_Hannum_1_1 <- EEAA_Hannum_1_1[!(row.names(EEAA_Hannum_1_1) %in% 
                                         c("2","11","12")),]
# Perform Shapiro test to check distribution
shapiro.test(x=c(EEAA_Hannum_1_1$baseline_EEAAHannum.ENmix,
                 EEAA_Hannum_1_1$baseline_EEAAHannumAge.Lisa)) 
# Perform paired T-test 
t.test(EEAA_Hannum_1_1$baseline_EEAAHannum.ENmix,
       EEAA_Hannum_1_1$baseline_EEAAHannumAge.Lisa,paired=T)

# Create dataframe with sample id and Visit 2 epigenetic age acceleration ENmix 
# values
Visit2_HannumENmix <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Sample_Name",
                                                       "EEAAHannum.ENmix")]
names(Visit2_HannumENmix)[2]<-"Visit2_EEAAHannum.ENmix"
# Create dataframe with sample id and Visit 2 epigenetic age acceleration Lisa
# values
Visit2_HannumLisa <- EEAAHannumAge[,c("Sample_Name","EEAAHannumAge")]
names(Visit2_HannumLisa)[2]<-"Visit2_EEAAHannumAge.Lisa"
# Merge into one dataframe with 3 columns
EEAA_Hannum_2_2 <- merge(Visit2_HannumENmix, Visit2_HannumLisa, 
                         by.x="Sample_Name",by.y="Sample_Name")
# Delete rows by name (delete those only on blockers: 10002, 10017, 10018)
EEAA_Hannum_2_2 <- EEAA_Hannum_2_2[!(row.names(EEAA_Hannum_2_2) %in% 
                                         c("2","11","12")),]
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Hannum_2_2$Visit2_EEAAHannum.ENmix,
                 EEAA_Hannum_2_2$Visit2_EEAAHannumAge.Lisa)) 
# Perform paired T-test
t.test(EEAA_Hannum_2_2$Visit2_EEAAHannum.ENmix,
       EEAA_Hannum_2_2$Visit2_EEAAHannumAge.Lisa,paired=T)


####### Check Distribution/Perform T-Tests for Visit 1/2 Age Acceleration #####
### All transgender youth on GAH:
# Horvath
# Create dataframe with sample id and baseline Horvath epigenetic age 
# acceleration ENmix values
baseline_EEAAHorvath <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                         "EEAAHorvath.ENmix")]
names(baseline_EEAAHorvath)[2]<-"baseline_EEAAHorvath"
# Create dataframe with sample id and Visit 2 Horvath epigenetic age 
# acceleration ENmix values
Visit2_EEAAHorvath <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                       "EEAAHorvath.ENmix")]
names(Visit2_EEAAHorvath)[2]<-"Visit2_EEAAHorvath"
# Merge into one dataframe with 3 columns
EEAA_Horvath_1_2 <- merge(baseline_EEAAHorvath, Visit2_EEAAHorvath, 
                          by.x="Subject.ID",by.y="Subject.ID")
# Delete rows by name (delete those only on blockers)
EEAA_Horvath_1_2 <- EEAA_Horvath_1_2[!(row.names(EEAA_Horvath_1_2) %in% 
                                           c("2","11","12")),]
saveRDS(EEAA_Horvath_1_2,file="G954_EPIC_Data_0217EEAA_Horvath_1_2_GAH.RDS")
write.csv(EEAA_Horvath_1_2,file="G954_EPIC_Data_0217EEAA_Horvath_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Horvath_1_2$baseline_EEAAHorvath,
                 EEAA_Horvath_1_2$Visit2_EEAAHorvath)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_Horvath_1_2$baseline_EEAAHorvath,
       EEAA_Horvath_1_2$Visit2_EEAAHorvath,paired=T)

## Hannum
# Create dataframe with sample id and baseline Hannum epigenetic age 
# acceleration ENmix values
baseline_EEAAHannum <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                        "EEAAHannum.ENmix")]
names(baseline_EEAAHannum)[2]<-"baseline_EEAAHannum"
# Create dataframe with sample id and Visit 2 Hannum epigenetic age 
# acceleration ENmix values
Visit2_EEAAHannum <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                      "EEAAHannum.ENmix")]
names(Visit2_EEAAHannum)[2]<-"Visit2_EEAAHannum"
# Merge into one dataframe with 3 columns
EEAA_Hannum_1_2 <- merge(baseline_EEAAHannum, Visit2_EEAAHannum, 
                         by.x="Subject.ID",by.y="Subject.ID")
# Delete rows by name (delete those only on blockers)
EEAA_Hannum_1_2 <- EEAA_Hannum_1_2[!(row.names(EEAA_Hannum_1_2) %in% 
                                         c("2","11","12")),]
saveRDS(EEAA_Hannum_1_2,file="G954_EPIC_Data_0217EEAA_Hannum_1_2_GAH.RDS")
write.csv(EEAA_Hannum_1_2,file="G954_EPIC_Data_0217EEAA_Hannum_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Hannum_1_2$baseline_EEAAHannum,
                 EEAA_Hannum_1_2$Visit2_EEAAHannum)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_Hannum_1_2$baseline_EEAAHannum,
       EEAA_Hannum_1_2$Visit2_EEAAHannum,paired=T)

## PhenoAge
# Create dataframe with sample id and baseline PhenoAge epigenetic age 
# acceleration ENmix values
baseline_EEAAPhenoAge <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                          "EEAAPhenoAge.ENmix")]
names(baseline_EEAAPhenoAge)[2]<-"baseline_EEAAPhenoAge"
# Create dataframe with sample id and Visit 2 PhenoAge epigenetic age 
# acceleration ENmix values
Visit2_EEAAPhenoAge <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                        "EEAAPhenoAge.ENmix")]
names(Visit2_EEAAPhenoAge)[2]<-"Visit2_EEAAPhenoAge"
# Merge into one dataframe with 3 columns
EEAA_PhenoAge_1_2 <- merge(baseline_EEAAPhenoAge, Visit2_EEAAPhenoAge, 
                           by.x="Subject.ID",by.y="Subject.ID")
# Delete rows by name (delete those only on blockers)
EEAA_PhenoAge_1_2 <- EEAA_PhenoAge_1_2[!(row.names(EEAA_PhenoAge_1_2) %in% 
                                             c("2","11","12")),]
saveRDS(EEAA_PhenoAge_1_2,file="G954_EPIC_Data_0217EEAA_PhenoAge_1_2_GAH.RDS")
write.csv(EEAA_PhenoAge_1_2,file="G954_EPIC_Data_0217EEAA_PhenoAge_1_2_GAH.csv")
# Perform Shaprio test to see distribution
shapiro.test(x=c(EEAA_PhenoAge_1_2$baseline_EEAAPhenoAge,
                 EEAA_PhenoAge_1_2$Visit2_EEAAPhenoAge)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_PhenoAge_1_2$baseline_EEAAPhenoAge,
       EEAA_PhenoAge_1_2$Visit2_EEAAPhenoAge,paired=T)

## GrimAge
# Create dataframe with sample id and baseline GrimAge epigenetic age 
# acceleration ENmix values
baseline_EEAAGrimAge <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                         "EEAAGrimAge")]
names(baseline_EEAAGrimAge)[2]<-"baseline_EEAAGrimAge"
# Create dataframe with sample id and Visit 2 GrimAge epigenetic age 
# acceleration ENmix values
Visit2_EEAAGrimAge <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                       "EEAAGrimAge")]
names(Visit2_EEAAGrimAge)[2]<-"Visit2_EEAAGrimAge"
# Merge into one dataframe with 3 columns
EEAA_GrimAge_1_2 <- merge(baseline_EEAAGrimAge, Visit2_EEAAGrimAge, 
                          by.x="Subject.ID",by.y="Subject.ID")
# Delete rows by name (delete those only on blockers)
EEAA_GrimAge_1_2 <- EEAA_GrimAge_1_2[!(row.names(EEAA_GrimAge_1_2) %in% 
                                           c("2","11","12")),]
saveRDS(EEAA_GrimAge_1_2,file="G954_EPIC_Data_0217EEAA_GrimAge_1_2_GAH.RDS")
write.csv(EEAA_GrimAge_1_2,file="G954_EPIC_Data_0217EEAA_GrimAge_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_GrimAge_1_2$baseline_EEAAGrimAge,
                 EEAA_GrimAge_1_2$Visit2_EEAAGrimAge)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_GrimAge_1_2$baseline_EEAAGrimAge,
       EEAA_GrimAge_1_2$Visit2_EEAAGrimAge,paired=T)


### All transgender youth on blockers:
## Horvath
# Create dataframe with sample id and baseline Horvath epigenetic age 
# acceleration ENmix values
baseline_EEAAHorvath <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                         "EEAAHorvath.ENmix")]
names(baseline_EEAAHorvath)[2]<-"baseline_EEAAHorvath"
# Create dataframe with sample id and Visit 2 Horvath epigenetic age 
# acceleration ENmix values
Visit2_EEAAHorvath <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                       "EEAAHorvath.ENmix")]
names(Visit2_EEAAHorvath)[2]<-"Visit2_EEAAHorvath"
# Merge into one dataframe with 3 columns
EEAA_Horvath_1_2 <- merge(baseline_EEAAHorvath, Visit2_EEAAHorvath, 
                          by.x="Subject.ID",by.y="Subject.ID")
# Select rows by name (select those only on blockers)
EEAA_Horvath_1_2 <- EEAA_Horvath_1_2[(row.names(EEAA_Horvath_1_2) %in% 
                                          c("2","11","12")),]
saveRDS(EEAA_Horvath_1_2,
        file="G954_EPIC_Data_0217EEAA_Horvath_1_2_Blockers.RDS")
write.csv(EEAA_Horvath_1_2,
          file="G954_EPIC_Data_0217EEAA_Horvath_1_2_Blockers.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Horvath_1_2$baseline_EEAAHorvath,
                 EEAA_Horvath_1_2$Visit2_EEAAHorvath)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_Horvath_1_2$baseline_EEAAHorvath,
       EEAA_Horvath_1_2$Visit2_EEAAHorvath,paired=T)

## Hannum
# Create dataframe with sample id and baseline Hannum epigenetic age 
# acceleration ENmix values
baseline_EEAAHannum <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                        "EEAAHannum.ENmix")]
names(baseline_EEAAHannum)[2]<-"baseline_EEAAHannum"
# Create dataframe with sample id and Visit 2 Hannum epigenetic age 
# acceleration ENmix values
Visit2_EEAAHannum <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                      "EEAAHannum.ENmix")]
names(Visit2_EEAAHannum)[2]<-"Visit2_EEAAHannum"
# Merge into one dataframe with 3 columns
EEAA_Hannum_1_2 <- merge(baseline_EEAAHannum, Visit2_EEAAHannum, 
                         by.x="Subject.ID",by.y="Subject.ID")
# Select rows by name (select those only on blockers)
EEAA_Hannum_1_2 <- EEAA_Hannum_1_2[(row.names(EEAA_Hannum_1_2) %in% 
                                        c("2","11","12")),]
saveRDS(EEAA_Hannum_1_2,file="G954_EPIC_Data_0217EEAA_Hannum_1_2_Blockers.RDS")
write.csv(EEAA_Hannum_1_2,file="G954_EPIC_Data_0217EEAA_Hannum_1_2_Blockers.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Hannum_1_2$baseline_EEAAHannum,
                 EEAA_Hannum_1_2$Visit2_EEAAHannum)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_Hannum_1_2$baseline_EEAAHannum,
       EEAA_Hannum_1_2$Visit2_EEAAHannum,paired=T)

## PhenoAge
# Create dataframe with sample id and baseline PhenoAge epigenetic age 
# acceleration ENmix values
baseline_EEAAPhenoAge <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                          "EEAAPhenoAge.ENmix")]
names(baseline_EEAAPhenoAge)[2]<-"baseline_EEAAPhenoAge"
# Create dataframe with sample id and Visit 2 PhenoAge epigenetic age 
# acceleration ENmix values
Visit2_EEAAPhenoAge <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                        "EEAAPhenoAge.ENmix")]
names(Visit2_EEAAPhenoAge)[2]<-"Visit2_EEAAPhenoAge"
# Merge into one dataframe with 3 columns
EEAA_PhenoAge_1_2 <- merge(baseline_EEAAPhenoAge, Visit2_EEAAPhenoAge, 
                           by.x="Subject.ID",by.y="Subject.ID")
# Select rows by name (select those only on blockers)
EEAA_PhenoAge_1_2 <- EEAA_PhenoAge_1_2[(row.names(EEAA_PhenoAge_1_2) %in% 
                                            c("2","11","12")),]
saveRDS(EEAA_PhenoAge_1_2,
        file="G954_EPIC_Data_0217EEAA_PhenoAge_1_2_Blockers.RDS")
write.csv(EEAA_PhenoAge_1_2,
          file="G954_EPIC_Data_0217EEAA_PhenoAge_1_2_Blockers.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_PhenoAge_1_2$baseline_EEAAPhenoAge,
                 EEAA_PhenoAge_1_2$Visit2_EEAAPhenoAge)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_PhenoAge_1_2$baseline_EEAAPhenoAge,
       EEAA_PhenoAge_1_2$Visit2_EEAAPhenoAge,paired=T)

## GrimAge
# Create dataframe with sample id and baseline GrimAge epigenetic age 
# acceleration ENmix values
baseline_EEAAGrimAge <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                         "EEAAGrimAge")]
names(baseline_EEAAGrimAge)[2]<-"baseline_EEAAGrimAge"
# Create dataframe with sample id and Visit 2 GrimAge epigenetic age 
# acceleration ENmix values
Visit2_EEAAGrimAge <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                       "EEAAGrimAge")]
names(Visit2_EEAAGrimAge)[2]<-"Visit2_EEAAGrimAge"
# Merge into one dataframe with 3 columns
EEAA_GrimAge_1_2 <- merge(baseline_EEAAGrimAge, Visit2_EEAAGrimAge, 
                          by.x="Subject.ID",by.y="Subject.ID")
# Select rows by name (select those only on blockers)
EEAA_GrimAge_1_2 <- EEAA_GrimAge_1_2[(row.names(EEAA_GrimAge_1_2) %in% 
                                          c("2","11","12")),]
saveRDS(EEAA_GrimAge_1_2,file="G954_EPIC_Data_0217EEAA_GrimAge_1_2_GAH.RDS")
write.csv(EEAA_GrimAge_1_2,file="G954_EPIC_Data_0217EEAA_GrimAge_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_GrimAge_1_2$baseline_EEAAGrimAge,
                 EEAA_GrimAge_1_2$Visit2_EEAAGrimAge)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_GrimAge_1_2$baseline_EEAAGrimAge,
       EEAA_GrimAge_1_2$Visit2_EEAAGrimAge,paired=T)


### All transmasculine youth on GAH:
## Horvath
# Create dataframe with sample id and baseline Horvath epigenetic age 
# acceleration ENmix values
baseline_EEAAHorvath <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                         "EEAAHorvath.ENmix")]
names(baseline_EEAAHorvath)[2]<-"baseline_EEAAHorvath"
# Create dataframe with sample id and Visit 2 Horvath epigenetic age 
# acceleration ENmix values
Visit2_EEAAHorvath <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                       "EEAAHorvath.ENmix")]
names(Visit2_EEAAHorvath)[2]<-"Visit2_EEAAHorvath"
# Merge into one dataframe with 3 columns
EEAA_Horvath_1_2 <- merge(baseline_EEAAHorvath, Visit2_EEAAHorvath, 
                          by.x="Subject.ID",by.y="Subject.ID")
# Delete rows by name (delete those only on blockers or assigned male at birth)
EEAA_Horvath_1_2 <- EEAA_Horvath_1_2[!(row.names(EEAA_Horvath_1_2) %in% 
                                           c("1","2","9","11","12","15")),]
saveRDS(EEAA_Horvath_1_2,file="G954_EPIC_Data_0217EEAA_Horvath_1_2_GAH.RDS")
write.csv(EEAA_Horvath_1_2,file="G954_EPIC_Data_0217EEAA_Horvath_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Horvath_1_2$baseline_EEAAHorvath,
                 EEAA_Horvath_1_2$Visit2_EEAAHorvath)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_Horvath_1_2$baseline_EEAAHorvath,
       EEAA_Horvath_1_2$Visit2_EEAAHorvath,paired=T)

## Hannum
# Create dataframe with sample id and baseline Hannum epigenetic age 
# acceleration ENmix values
baseline_EEAAHannum <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                        "EEAAHannum.ENmix")]
names(baseline_EEAAHannum)[2]<-"baseline_EEAAHannum"
# Create dataframe with sample id and Visit 2 Hannum epigenetic age 
# acceleration ENmix values
Visit2_EEAAHannum <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                      "EEAAHannum.ENmix")]
names(Visit2_EEAAHannum)[2]<-"Visit2_EEAAHannum"
# Merge into one dataframe with 3 columns
EEAA_Hannum_1_2 <- merge(baseline_EEAAHannum, Visit2_EEAAHannum, 
                         by.x="Subject.ID",by.y="Subject.ID")
# Delete rows by name (delete those only on blockers or assigned male at birth)
EEAA_Hannum_1_2 <- EEAA_Hannum_1_2[!(row.names(EEAA_Hannum_1_2) %in% 
                                         c("1","2","9","11","12","15")),]
saveRDS(EEAA_Hannum_1_2,file="G954_EPIC_Data_0217EEAA_Hannum_1_2_GAH.RDS")
write.csv(EEAA_Hannum_1_2,file="G954_EPIC_Data_0217EEAA_Hannum_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Hannum_1_2$baseline_EEAAHannum,
                 EEAA_Hannum_1_2$Visit2_EEAAHannum)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_Hannum_1_2$baseline_EEAAHannum,
       EEAA_Hannum_1_2$Visit2_EEAAHannum,paired=T)

## PhenoAge
# Create dataframe with sample id and baseline PhenoAge epigenetic age 
# acceleration ENmix values
baseline_EEAAPhenoAge <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                          "EEAAPhenoAge.ENmix")]
names(baseline_EEAAPhenoAge)[2]<-"baseline_EEAAPhenoAge"
# Create dataframe with sample id and Visit 2 PhenoAge epigenetic age 
# acceleration ENmix values
Visit2_EEAAPhenoAge <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                        "EEAAPhenoAge.ENmix")]
names(Visit2_EEAAPhenoAge)[2]<-"Visit2_EEAAPhenoAge"
# Merge into one dataframe with 3 columns
EEAA_PhenoAge_1_2 <- merge(baseline_EEAAPhenoAge, Visit2_EEAAPhenoAge, 
                           by.x="Subject.ID",by.y="Subject.ID")
# Delete rows by name (delete those only on blockers or assigned male at birth)
EEAA_PhenoAge_1_2 <- EEAA_PhenoAge_1_2[!(row.names(EEAA_PhenoAge_1_2) %in% 
                                             c("1","2","9","11","12","15")),]
saveRDS(EEAA_PhenoAge_1_2,file="G954_EPIC_Data_0217EEAA_PhenoAge_1_2_GAH.RDS")
write.csv(EEAA_PhenoAge_1_2,file="G954_EPIC_Data_0217EEAA_PhenoAge_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_PhenoAge_1_2$baseline_EEAAPhenoAge,
                 EEAA_PhenoAge_1_2$Visit2_EEAAPhenoAge)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_PhenoAge_1_2$baseline_EEAAPhenoAge,
       EEAA_PhenoAge_1_2$Visit2_EEAAPhenoAge,paired=T)

## GrimAge
# Create dataframe with sample id and baseline GrimAge epigenetic age 
# acceleration ENmix values
baseline_EEAAGrimAge <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                         "EEAAGrimAge")]
names(baseline_EEAAGrimAge)[2]<-"baseline_EEAAGrimAge"
# Create dataframe with sample id and Visit 2 GrimAge epigenetic age 
# acceleration ENmix values
Visit2_EEAAGrimAge <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                       "EEAAGrimAge")]
names(Visit2_EEAAGrimAge)[2]<-"Visit2_EEAAGrimAge"
# Merge into one dataframe with 3 columns
EEAA_GrimAge_1_2 <- merge(baseline_EEAAGrimAge, Visit2_EEAAGrimAge, 
                          by.x="Subject.ID",by.y="Subject.ID")
# Delete rows by name (delete those only on blockers or assigned male at birth)
EEAA_GrimAge_1_2 <- EEAA_GrimAge_1_2[!(row.names(EEAA_GrimAge_1_2) %in% 
                                           c("1","2","9","11","12","15")),]
saveRDS(EEAA_GrimAge_1_2,file="G954_EPIC_Data_0217EEAA_GrimAge_1_2_GAH.RDS")
write.csv(EEAA_GrimAge_1_2,file="G954_EPIC_Data_0217EEAA_GrimAge_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_GrimAge_1_2$baseline_EEAAGrimAge,
                 EEAA_GrimAge_1_2$Visit2_EEAAGrimAge)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_GrimAge_1_2$baseline_EEAAGrimAge,
       EEAA_GrimAge_1_2$Visit2_EEAAGrimAge,paired=T)


### All transfeminine youth on GAH:
## Horvath
# Create dataframe with sample id and baseline Horvath epigenetic age 
# acceleration ENmix values
baseline_EEAAHorvath <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                         "EEAAHorvath.ENmix")]
names(baseline_EEAAHorvath)[2]<-"baseline_EEAAHorvath"
# Create dataframe with sample id and Visit 2 Horvath epigenetic age 
# acceleration ENmix values
Visit2_EEAAHorvath <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                       "EEAAHorvath.ENmix")]
names(Visit2_EEAAHorvath)[2]<-"Visit2_EEAAHorvath"
# Merge into one dataframe with 3 columns
EEAA_Horvath_1_2 <- merge(baseline_EEAAHorvath, Visit2_EEAAHorvath, 
                          by.x="Subject.ID",by.y="Subject.ID")
# Select rows by name (select those only assigned male at birth)
EEAA_Horvath_1_2 <- EEAA_Horvath_1_2[(row.names(EEAA_Horvath_1_2) %in% 
                                          c("1","9","15")),]
saveRDS(EEAA_Horvath_1_2,
        file="G954_EPIC_Data_0217EEAA_Horvath_1_2_Blockers.RDS")
write.csv(EEAA_Horvath_1_2,
          file="G954_EPIC_Data_0217EEAA_Horvath_1_2_Blockers.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Horvath_1_2$baseline_EEAAHorvath,
                 EEAA_Horvath_1_2$Visit2_EEAAHorvath)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_Horvath_1_2$baseline_EEAAHorvath,
       EEAA_Horvath_1_2$Visit2_EEAAHorvath,paired=T)

## Hannum
# Create dataframe with sample id and baseline Hannum epigenetic age 
# acceleration ENmix values
baseline_EEAAHannum <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                        "EEAAHannum.ENmix")]
names(baseline_EEAAHannum)[2]<-"baseline_EEAAHannum"
# Create dataframe with sample id and Visit 2 Hannum epigenetic age 
# acceleration ENmix values
Visit2_EEAAHannum <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                      "EEAAHannum.ENmix")]
names(Visit2_EEAAHannum)[2]<-"Visit2_EEAAHannum"
# Merge into one dataframe with 3 columns
EEAA_Hannum_1_2 <- merge(baseline_EEAAHannum, Visit2_EEAAHannum, 
                         by.x="Subject.ID",by.y="Subject.ID")
# Select rows by name (select those only assigned male at birth)
EEAA_Hannum_1_2 <- EEAA_Hannum_1_2[(row.names(EEAA_Hannum_1_2) %in% 
                                        c("1","9","15")),]
saveRDS(EEAA_Hannum_1_2,file="G954_EPIC_Data_0217EEAA_Hannum_1_2_Blockers.RDS")
write.csv(EEAA_Hannum_1_2,file="G954_EPIC_Data_0217EEAA_Hannum_1_2_Blockers.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_Hannum_1_2$baseline_EEAAHannum,
                 EEAA_Hannum_1_2$Visit2_EEAAHannum)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_Hannum_1_2$baseline_EEAAHannum,
       EEAA_Hannum_1_2$Visit2_EEAAHannum,paired=T)

## PhenoAge
# Create dataframe with sample id and baseline PhenoAge epigenetic age 
# acceleration ENmix values
baseline_EEAAPhenoAge <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                          "EEAAPhenoAge.ENmix")]
names(baseline_EEAAPhenoAge)[2]<-"baseline_EEAAPhenoAge"
# Create dataframe with sample id and Visit 2 PhenoAge epigenetic age 
# acceleration ENmix values
Visit2_EEAAPhenoAge <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                        "EEAAPhenoAge.ENmix")]
names(Visit2_EEAAPhenoAge)[2]<-"Visit2_EEAAPhenoAge"
# Merge into one dataframe with 3 columns
EEAA_PhenoAge_1_2 <- merge(baseline_EEAAPhenoAge, Visit2_EEAAPhenoAge, 
                           by.x="Subject.ID",by.y="Subject.ID")
# Select rows by name (select those only assigned male at birth)
EEAA_PhenoAge_1_2 <- EEAA_PhenoAge_1_2[(row.names(EEAA_PhenoAge_1_2) %in% 
                                            c("1","9","15")),]
saveRDS(EEAA_PhenoAge_1_2,
        file="G954_EPIC_Data_0217EEAA_PhenoAge_1_2_Blockers.RDS")
write.csv(EEAA_PhenoAge_1_2,
          file="G954_EPIC_Data_0217EEAA_PhenoAge_1_2_Blockers.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_PhenoAge_1_2$baseline_EEAAPhenoAge,
                 EEAA_PhenoAge_1_2$Visit2_EEAAPhenoAge)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_PhenoAge_1_2$baseline_EEAAPhenoAge,
       EEAA_PhenoAge_1_2$Visit2_EEAAPhenoAge,paired=T)

## GrimAge
# Create dataframe with sample id and baseline GrimAge epigenetic age 
# acceleration ENmix values
baseline_EEAAGrimAge <- EEAA.ENmix[EEAA.ENmix$Visit==1,c("Subject.ID",
                                                         "EEAAGrimAge")]
names(baseline_EEAAGrimAge)[2]<-"baseline_EEAAGrimAge"
# Create dataframe with sample id and Visit 2 GrimAge epigenetic age 
# acceleration ENmix values
Visit2_EEAAGrimAge <- EEAA.ENmix[EEAA.ENmix$Visit==2,c("Subject.ID",
                                                       "EEAAGrimAge")]
names(Visit2_EEAAGrimAge)[2]<-"Visit2_EEAAGrimAge"
# Merge into one dataframe with 3 columns
EEAA_GrimAge_1_2 <- merge(baseline_EEAAGrimAge, Visit2_EEAAGrimAge, 
                          by.x="Subject.ID",by.y="Subject.ID")
# Select rows by name (select those only assigned male at birth)
EEAA_GrimAge_1_2 <- EEAA_GrimAge_1_2[(row.names(EEAA_GrimAge_1_2) %in% 
                                          c("1","9","15")),]
saveRDS(EEAA_GrimAge_1_2,file="G954_EPIC_Data_0217EEAA_GrimAge_1_2_GAH.RDS")
write.csv(EEAA_GrimAge_1_2,file="G954_EPIC_Data_0217EEAA_GrimAge_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(EEAA_GrimAge_1_2$baseline_EEAAGrimAge,
                 EEAA_GrimAge_1_2$Visit2_EEAAGrimAge)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(EEAA_GrimAge_1_2$baseline_EEAAGrimAge,
       EEAA_GrimAge_1_2$Visit2_EEAAGrimAge,paired=T)