# Transgender Congruence Scale

# Load packages
library(dplyr)
library(piggyback)

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data")

# Load data
TCSdata <- read.csv("TCSdata.csv")

# Replace skipped questions (response '9') with NA
TCSdata %>% 
    replace(., .==9,NA) ->
    TCSdata

# Reverse-score columns
TCSdata$nrepresent <- recode(TCSdata$nrepresent,
                                              "1"=5, "2"=4, "3"=3, "4"=2, "5"=1)
TCSdata$nreflect <- recode(TCSdata$nreflect,
                           "1"=5, "2"=4, "3"=3, "4"=2, "5"=1)
TCSdata$nproud <- recode(TCSdata$nproud,
                         "1"=5, "2"=4, "3"=3, "4"=2, "5"=1)

# Make new column with row means of TCS scores
TCSdata$TCS <- rowMeans(subset(TCSdata,
                                     select=c("outward",
                                              "unity",
                                              "adequate",
                                              "comfort",
                                              "physical",
                                              "nrepresent",
                                              "happyapp",
                                              "nreflect",
                                              "consistent",
                                              "nproud",
                                              "happygid",
                                              "acceptgid")),na.rm=TRUE)
# Make new column with row means of Appearance Congruence Subscale scores
TCSdata$ACS <- rowMeans(subset(TCSdata, 
                                     select=c("outward",
                                              "unity",
                                              "adequate",
                                              "comfort",
                                              "physical",
                                              "nrepresent",
                                              "happyapp",
                                              "nreflect",
                                              "consistent")),na.rm=TRUE)
# Make new column with row means of Gender Identity Acceptance Subscale scores
TCSdata$GIAS <- rowMeans(subset(TCSdata, 
                                      select=c("nproud",
                                               "happygid",
                                               "acceptgid")),na.rm=TRUE)


# Move score columns to front
TCSdata %>% relocate(Subject.ID, Visit, TCS, ACS, GIAS) ->
    TCSdata

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
write.csv(TCSdata, file="TCSdata.csv")

#### T-Tests ####

## All transgender youth on GAH
# Create dataframe with sample id and baseline TCS values
baseline_TCS <- TCSdata[TCSdata$Visit==1,c("Subject.ID","TCS")]
names(baseline_TCS)[2]<-"baseline_TCS"
# Create dataframe with sample id and Visit 2 TCS values
Visit2_TCS <- TCSdata[TCSdata$Visit==2,c("Subject.ID","TCS")]
names(Visit2_TCS)[2]<-"Visit2_TCS"
# Merge into one dataframe with 3 columns
TCS_1_2 <- merge(baseline_TCS, Visit2_TCS, by.x="Subject.ID",by.y="Subject.ID",
                 all.x=F)
# Delete rows by name (delete those only on blockers)
TCS_1_2 <- TCS_1_2[!(row.names(TCS_1_2) %in% c("2","11","12")),]
saveRDS(TCS_1_2,file="TCS_1_2_GAH.RDS")
write.csv(TCS_1_2,file="TCS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(TCS_1_2$baseline_TCS,
                 TCS_1_2$Visit2_TCS)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(TCS_1_2$baseline_TCS,
       TCS_1_2$Visit2_TCS,paired=T)

# Create dataframe with sample id and baseline ACS values
baseline_ACS <- TCSdata[TCSdata$Visit==1,c("Subject.ID","ACS")]
names(baseline_ACS)[2]<-"baseline_ACS"
# Create dataframe with sample id and Visit 2 ACS values
Visit2_ACS <- TCSdata[TCSdata$Visit==2,c("Subject.ID","ACS")]
names(Visit2_ACS)[2]<-"Visit2_ACS"
# Merge into one dataframe with 3 columns
ACS_1_2 <- merge(baseline_ACS, Visit2_ACS, by.x="Subject.ID",by.y="Subject.ID",
                 all.x=F)
# Delete rows by name (delete those only on blockers)
ACS_1_2 <- ACS_1_2[!(row.names(ACS_1_2) %in% c("2","11","12")),]
saveRDS(ACS_1_2,file="ACS_1_2_GAH.RDS")
write.csv(ACS_1_2,file="ACS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(ACS_1_2$baseline_ACS,
                 ACS_1_2$Visit2_ACS)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(ACS_1_2$baseline_ACS,
       ACS_1_2$Visit2_ACS,paired=T)

# Create dataframe with sample id and baseline GIAS  values
baseline_GIAS <- TCSdata[TCSdata$Visit==1,c("Subject.ID","GIAS")]
names(baseline_GIAS)[2]<-"baseline_GIAS"
# Create dataframe with sample id and Visit 2 GIAS values
Visit2_GIAS <- TCSdata[TCSdata$Visit==2,c("Subject.ID","GIAS")]
names(Visit2_GIAS)[2]<-"Visit2_GIAS"
# Merge into one dataframe with 3 columns
GIAS_1_2 <- merge(baseline_GIAS, Visit2_GIAS, by.x="Subject.ID",
                  by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
GIAS_1_2 <- GIAS_1_2[!(row.names(GIAS_1_2) %in% c("2","11","12")),]
saveRDS(GIAS_1_2,file="GIAS_1_2_GAH.RDS")
write.csv(GIAS_1_2,file="GIAS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(GIAS_1_2$baseline_GIAS,
                 GIAS_1_2$Visit2_GIAS)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(GIAS_1_2$baseline_GIAS,
       GIAS_1_2$Visit2_GIAS,paired=T)


## Transmasculine youth on GAH
# Create dataframe with sample id and baseline TCS values
baseline_TCS <- TCSdata[TCSdata$Visit==1,c("Subject.ID","TCS")]
names(baseline_TCS)[2]<-"baseline_TCS"
# Create dataframe with sample id and Visit 2 TCS values
Visit2_TCS <- TCSdata[TCSdata$Visit==2,c("Subject.ID","TCS")]
names(Visit2_TCS)[2]<-"Visit2_TCS"
# Merge into one dataframe with 3 columns
TCS_1_2 <- merge(baseline_TCS, Visit2_TCS, by.x="Subject.ID",by.y="Subject.ID",
                 all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
TCS_1_2 <- TCS_1_2[!(row.names(TCS_1_2) %in% c("1","2","9","11","12","16")),]
saveRDS(TCS_1_2,file="TCS_1_2_GAH.RDS")
write.csv(TCS_1_2,file="TCS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(TCS_1_2$baseline_TCS,
                 TCS_1_2$Visit2_TCS)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(TCS_1_2$baseline_TCS,
       TCS_1_2$Visit2_TCS,paired=T)


# Create dataframe with sample id and baseline ACS values
baseline_ACS <- TCSdata[TCSdata$Visit==1,c("Subject.ID","ACS")]
names(baseline_ACS)[2]<-"baseline_ACS"
# Create dataframe with sample id and Visit 2 ACS values
Visit2_ACS <- TCSdata[TCSdata$Visit==2,c("Subject.ID","ACS")]
names(Visit2_ACS)[2]<-"Visit2_ACS"
# Merge into one dataframe with 3 columns
ACS_1_2 <- merge(baseline_ACS, Visit2_ACS, by.x="Subject.ID",by.y="Subject.ID",
                 all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
ACS_1_2 <- ACS_1_2[!(row.names(ACS_1_2) %in% c("1","2","9","11","12","16")),]
saveRDS(ACS_1_2,file="ACS_1_2_GAH.RDS")
write.csv(ACS_1_2,file="ACS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(ACS_1_2$baseline_ACS,
                 ACS_1_2$Visit2_ACS)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(ACS_1_2$baseline_ACS,
       ACS_1_2$Visit2_ACS,paired=T)

# Create dataframe with sample id and baseline GIAS values
baseline_GIAS <- TCSdata[TCSdata$Visit==1,c("Subject.ID","GIAS")]
names(baseline_GIAS)[2]<-"baseline_GIAS"
# Create dataframe with sample id and Visit 2 GIAS values
Visit2_GIAS <- TCSdata[TCSdata$Visit==2,c("Subject.ID","GIAS")]
names(Visit2_GIAS)[2]<-"Visit2_GIAS"
# Merge into one dataframe with 3 columns
GIAS_1_2 <- merge(baseline_GIAS, Visit2_GIAS, by.x="Subject.ID",
                  by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
GIAS_1_2 <- GIAS_1_2[!(row.names(GIAS_1_2) %in% c("1","2","9","11","12","16")),]
saveRDS(GIAS_1_2,file="GIAS_1_2_GAH.RDS")
write.csv(GIAS_1_2,file="GIAS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(GIAS_1_2$baseline_GIAS,
                 GIAS_1_2$Visit2_GIAS)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(GIAS_1_2$baseline_GIAS,
       GIAS_1_2$Visit2_GIAS,paired=T)