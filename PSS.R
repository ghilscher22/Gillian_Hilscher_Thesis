# 10-Item Perceived Stress Scale

# Load packages
library(dplyr)
library(piggyback)

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data")

# Load data
PSSdata <- read.csv("PSSdata.csv")

# Replace skipped questions (response '15') with NA
PSSdata %>% 
    replace(., .==15,NA) ->
    PSSdata

# Reverse-score columns 
PSSdata$personalprob <- recode(PSSdata$personalprob,
                                      "0"=4, "1"=3, "2"=2, "3"=1, "4"=0)
PSSdata$myway <- recode(PSSdata$myway,
                                      "0"=4, "1"=3, "2"=2, "3"=1, "4"=0)
PSSdata$irritation <- recode(PSSdata$irritation,
                             "0"=4, "1"=3, "2"=2, "3"=1, "4"=0)
PSSdata$on_top <- recode(PSSdata$on_top,
                         "0"=4, "1"=3, "2"=2, "3"=1, "4"=0)

# Make new columns with row means
PSSdata$PSS <- rowMeans(subset(PSSdata, 
                               select=c("unexpected",
                                        "control",
                                        "nervousstress",
                                        "personalprob",
                                        "myway",
                                        "cope",
                                        "irritation",
                                        "on_top",
                                        "nocontrol_anger",
                                        "noovercome")),na.rm=TRUE)

# Move score columns to front
PSSdata %>% relocate(Subject.ID, Visit, PSS) -> 
    PSSdata

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
write.csv(PSSdata,file="PSSdata.csv")


#### T-Tests ####

## All transgender youth on GAH
# Create dataframe with sample id and baseline PSS values
baseline_PSS <- PSSdata[PSSdata$Visit==1,c("Subject.ID","PSS")]
names(baseline_PSS)[2]<-"baseline_PSS"
# Create dataframe with sample id and Visit 2 PSS values
Visit2_PSS <- PSSdata[PSSdata$Visit==2,c("Subject.ID","PSS")]
names(Visit2_PSS)[2]<-"Visit2_PSS"
# Merge into one dataframe with 3 columns
PSS_1_2 <- merge(baseline_PSS, Visit2_PSS, by.x="Subject.ID",by.y="Subject.ID",
                 all.x=F)
# Delete rows by name (delete those only on blockers)
PSS_1_2 <- PSS_1_2[!(row.names(PSS_1_2) %in% c("2","11","12")),]
saveRDS(PSS_1_2,file="PSS_1_2_GAH.RDS")
write.csv(PSS_1_2,file="PSS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(PSS_1_2$baseline_PSS,
                 PSS_1_2$Visit2_PSS))
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(PSS_1_2$baseline_PSS,
       PSS_1_2$Visit2_PSS,paired=T)


## Transmasculine youth on GAH
# Create dataframe with sample id and baseline PSS values
baseline_PSS <- PSSdata[PSSdata$Visit==1,c("Subject.ID","PSS")]
names(baseline_PSS)[2]<-"baseline_PSS"
# Create dataframe with sample id and Visit 2 PSS values
Visit2_PSS <- PSSdata[PSSdata$Visit==2,c("Subject.ID","PSS")]
names(Visit2_PSS)[2]<-"Visit2_PSS"
# Merge into one dataframe with 3 columns
PSS_1_2 <- merge(baseline_PSS, Visit2_PSS, by.x="Subject.ID",by.y="Subject.ID",
                 all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
PSS_1_2 <- PSS_1_2[!(row.names(PSS_1_2) %in% c("1","2","9","11","12","16")),]
saveRDS(PSS_1_2,file="PSS_1_2_GAH.RDS")
write.csv(PSS_1_2,file="PSS_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(PSS_1_2$baseline_PSS,
                 PSS_1_2$Visit2_PSS))
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(PSS_1_2$baseline_PSS,
       PSS_1_2$Visit2_PSS,paired=T) 