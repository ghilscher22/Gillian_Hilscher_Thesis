# Objectified Body Consciousness Scale: Surveillance Subscale

# Load packages
library(dplyr)
library(piggyback)

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data")

# Load data
OBCSdata <- read.csv("OBCSdata.csv")

# Replace responses of "skip" (9) or "NA" (8) with NA
OBCSdata %>% 
    replace(., .==8 | .==9,NA) ->
    OBCSdata

# Make new columns with row means
OBCSdata$OBCS <- rowMeans(subset(OBCSdata,
                                 select=c("subject_often",
                                          "subject_many",
                                          "subject_good",
                                          "subject_other")),na.rm=TRUE)

# Move score columns to front
OBCSdata %>% relocate(Subject.ID, Visit, OBCS) ->
    OBCSdata

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
write.csv(OBCSdata,file="OBCSdata.csv")

#### T-Tests ####

## All transgender youth on GAH
# Create dataframe with sample id and baseline OBCS values
baseline_OBCS <- OBCSdata[OBCSdata$Visit==1,c("Subject.ID","OBCS")]
names(baseline_OBCS)[2]<-"baseline_OBCS"
# Create dataframe with sample id and Visit 2 OBCS values
Visit2_OBCS <- OBCSdata[OBCSdata$Visit==2,c("Subject.ID","OBCS")]
names(Visit2_OBCS)[2]<-"Visit2_OBCS"
# Merge into one dataframe with 3 columns
OBCS_1_2 <- merge(baseline_OBCS, Visit2_OBCS, by.x="Subject.ID",
                  by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
OBCS_1_2 <- OBCS_1_2[!(row.names(OBCS_1_2) %in% c("2","11","12")),]
saveRDS(OBCS_1_2,file="OBCS_1_2_GAH.RDS")
write.csv(OBCS_1_2,file="OBCS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(OBCS_1_2$baseline_OBCS,
                 OBCS_1_2$Visit2_OBCS)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(OBCS_1_2$baseline_OBCS,
       OBCS_1_2$Visit2_OBCS,paired=T)


## Transmasculine youth on GAH
# Create dataframe with sample id and baseline OBCS values
baseline_OBCS <- OBCSdata[OBCSdata$Visit==1,c("Subject.ID","OBCS")]
names(baseline_OBCS)[2]<-"baseline_OBCS"
# Create dataframe with sample id and Visit 2 OBCS values
Visit2_OBCS <- OBCSdata[OBCSdata$Visit==2,c("Subject.ID","OBCS")]
names(Visit2_OBCS)[2]<-"Visit2_OBCS"
# Merge into one dataframe with 3 columns
OBCS_1_2 <- merge(baseline_OBCS, Visit2_OBCS, by.x="Subject.ID",
                  by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
OBCS_1_2 <- OBCS_1_2[!(row.names(OBCS_1_2) %in% c("1","2","9","11","12","16")),]
saveRDS(OBCS_1_2,file="OBCS_1_2_GAH.RDS")
write.csv(OBCS_1_2,file="OBCS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(OBCS_1_2$baseline_OBCS,
                 OBCS_1_2$Visit2_OBCS)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(OBCS_1_2$baseline_OBCS,
       OBCS_1_2$Visit2_OBCS,paired=T)