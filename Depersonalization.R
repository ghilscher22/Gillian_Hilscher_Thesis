# Depersonalization

# Load packages
library(dplyr)
library(piggyback)

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data")

# Load data
Depersdata <- read.csv("Depersonalizationdata.csv")

# Replace skipped questions (response '8' or '9') with NA
Depersdata %>% 
    replace(., .==8 | .==9,NA) ->
    Depersdata

# Make new columns with row means
Depersdata$Depersonalization <- rowMeans(subset(Depersdata, 
                                              select=c("subject_appchange",
                                              "subject_nfood",
                                              "subject_ashamedphys",
                                              "subject_mirror",
                                              "subject_detached",
                                              "subject_nbelong")),na.rm=TRUE)

# Move score columns to front
Depersdata %>% relocate(Subject.ID, Visit, Depersonalization) ->
    Depersdata

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
write.csv(Depersdata,file="Depersdata.csv")

#### T-Tests ####

## All transgender youth on GAH
# Create dataframe with sample id and baseline Depersonalization values
baseline_Depers <- Depersdata[Depersdata$Visit==1,
                              c("Subject.ID","Depersonalization")]
names(baseline_Depers)[2]<-"baseline_Depers"
# Create dataframe with sample id and Visit 2 Depersonalization values
Visit2_Depers <- Depersdata[Depersdata$Visit==2,
                            c("Subject.ID","Depersonalization")]
names(Visit2_Depers)[2]<-"Visit2_Depers"
# Merge into one dataframe with 3 columns
Depers_1_2 <- merge(baseline_Depers, Visit2_Depers, by.x="Subject.ID",
                    by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
Depers_1_2 <- Depers_1_2[!(row.names(Depers_1_2) %in% c("2","11","12")),]
saveRDS(Depers_1_2,file="Depers_1_2_GAH.RDS")
write.csv(Depers_1_2,file="Depers_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(Depers_1_2$baseline_Depers,
                 Depers_1_2$Visit2_Depers)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(Depers_1_2$baseline_Depers,
       Depers_1_2$Visit2_Depers,paired=T)


## Transmasculine youth on GAH
# Create dataframe with sample id and baseline Depersonalization values
baseline_Depers <- Depersdata[Depersdata$Visit==1,
                              c("Subject.ID","Depersonalization")]
names(baseline_Depers)[2]<-"baseline_Depers"
# Create dataframe with sample id and Visit 2 Depersonalization values
Visit2_Depers <- Depersdata[Depersdata$Visit==2,
                            c("Subject.ID","Depersonalization")]
names(Visit2_Depers)[2]<-"Visit2_Depers"
# Merge into one dataframe with 3 columns
Depers_1_2 <- merge(baseline_Depers, Visit2_Depers, by.x="Subject.ID",
                    by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
Depers_1_2 <- Depers_1_2[!(row.names(Depers_1_2) %in% 
                               c("1","2","9","11","12","16")),]
saveRDS(Depers_1_2,file="Depers_1_2_GAH.RDS")
write.csv(Depers_1_2,file="Depers_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(Depers_1_2$baseline_Depers,
                 Depers_1_2$Visit2_Depers)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(Depers_1_2$baseline_Depers,
       Depers_1_2$Visit2_Depers,paired=T)