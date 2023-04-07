# Rosenberg Self-Esteem Scale

# Load packages
library(dplyr)
library(piggyback)

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data")

# Load data
Selfesteemdata <- read.csv("Selfesteemdata.csv")

# Replace skipped questions (response '9') with NA
Selfesteemdata %>% 
    replace(., .==9,NA) ->
    Selfesteemdata

# Reverse-score columns 
Selfesteemdata$subject_nogood <- recode(Selfesteemdata$subject_nogood,
                                      "1"=4, "2"=3, "3"=2, "4"=1)
Selfesteemdata$subject_muchpride <- recode(Selfesteemdata$subject_muchpride,
                                           "1"=4, "2"=3, "3"=2, "4"=1)
Selfesteemdata$subject_useless <- recode(Selfesteemdata$subject_useless,
                                         "1"=4, "2"=3, "3"=2, "4"=1)
Selfesteemdata$subject_wishrespect <- recode(Selfesteemdata$subject_wishrespect,
                                             "1"=4, "2"=3, "3"=2, "4"=1)
Selfesteemdata$subject_failure <- recode(Selfesteemdata$subject_failure,
                                         "1"=4, "2"=3, "3"=2, "4"=1)

# Make new columns with row means
Selfesteemdata$Self.Esteem <- rowMeans(subset(Selfesteemdata, 
                                            select=c("subject_satisfied",
                                            "subject_nogood",
                                            "subject_gqual",
                                            "subject_dowell",
                                            "subject_muchpride",
                                            "subject_useless",
                                            "subject_worth",
                                            "subject_wishrespect",
                                            "subject_failure",
                                            "subject_positive")),na.rm=TRUE)

# mMve score columns to front
Selfesteemdata %>% 
    relocate(Subject.ID, Visit, Self.Esteem) ->
    Selfesteemdata

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data/03-24-23")
write.csv(Selfesteemdata,file="Selfesteemdata.csv")


#### T-Tests ####

## All transgender youth on GAH
# Create dataframe with sample id and baseline self-esteem values
baseline_Selfesteem <- Selfesteemdata[Selfesteemdata$Visit==1,
                                      c("Subject.ID","Self.Esteem")]
names(baseline_Selfesteem)[2]<-"baseline_Selfesteem"
# Create dataframe with sample id and Visit 2 self-esteem values
Visit2_Selfesteem <- Selfesteemdata[Selfesteemdata$Visit==2,
                                    c("Subject.ID","Self.Esteem")]
names(Visit2_Selfesteem)[2]<-"Visit2_Selfesteem"
# Merge into one dataframe with 3 columns
Selfesteem_1_2 <- merge(baseline_Selfesteem, Visit2_Selfesteem, 
                        by.x="Subject.ID",by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
Selfesteem_1_2 <- Selfesteem_1_2[!(row.names(Selfesteem_1_2) %in% 
                                       c("2","11","12")),]
saveRDS(Selfesteem_1_2,file="Selfesteem_1_2_GAH.RDS")
write.csv(Selfesteem_1_2,file="Selfesteem_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(Selfesteem_1_2$baseline_Selfesteem,
                 Selfesteem_1_2$Visit2_Selfesteem)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(Selfesteem_1_2$baseline_Selfesteem,
       Selfesteem_1_2$Visit2_Selfesteem,paired=T)


## Transmasculine youth on GAH
# Create dataframe with sample id and baseline self-esteem  values
baseline_Selfesteem <- Selfesteemdata[Selfesteemdata$Visit==1,
                                      c("Subject.ID","Self.Esteem")]
names(baseline_Selfesteem)[2]<-"baseline_Selfesteem"
# Create dataframe with sample id and Visit 2 self-esteem values
Visit2_Selfesteem <- Selfesteemdata[Selfesteemdata$Visit==2,
                                    c("Subject.ID","Self.Esteem")]
names(Visit2_Selfesteem)[2]<-"Visit2_Selfesteem"
# Merge into one dataframe with 3 columns
Selfesteem_1_2 <- merge(baseline_Selfesteem, Visit2_Selfesteem, 
                        by.x="Subject.ID",by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
Selfesteem_1_2 <- Selfesteem_1_2[!(row.names(Selfesteem_1_2) %in% 
                                       c("1","2","9","11","12","16")),]
saveRDS(Selfesteem_1_2,file="Selfesteem_1_2_GAH.RDS")
write.csv(Selfesteem_1_2,file="Selfesteem_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(Selfesteem_1_2$baseline_Selfesteem,
                 Selfesteem_1_2$Visit2_Selfesteem)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(Selfesteem_1_2$baseline_Selfesteem,
       Selfesteem_1_2$Visit2_Selfesteem,paired=T)