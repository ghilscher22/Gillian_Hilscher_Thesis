# PROMIS Anxiety and Depression

# Load packages
library(dplyr)
library(piggyback)
library(tidyverse)

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data")

# Load data
PROMISdata <- read.csv("PROMISdata.csv")

# Rescore all pediatric responses from 0-4 to 1-5
PROMISdata %>% 
    mutate(promawful = ifelse(!is.na(promawful) & 
                                  promawful == 15, NA,
                              ifelse(!is.na(promawful), 
                                     promawful + 1, NA)),
           promnervous = ifelse(!is.na(promnervous) & 
                                    promnervous == 15, NA,
                                      ifelse(!is.na(promnervous), 
                                             promnervous + 1, NA)),
           promscared = ifelse(!is.na(promscared) & 
                                   promscared == 15, NA,
                                      ifelse(!is.na(promscared), 
                                             promscared + 1, NA)),
           promworried = ifelse(!is.na(promworried) & 
                                    promworried == 15, NA,
                                      ifelse(!is.na(promworried), 
                                             promworried + 1, NA)),
           promworryhome = ifelse(!is.na(promworryhome) & 
                                      promworryhome == 15, NA,
                                      ifelse(!is.na(promworryhome), 
                                             promworryhome + 1, NA)),
           promeasyscare = ifelse(!is.na(promeasyscare) & 
                                      promeasyscare == 15, NA,
                                      ifelse(!is.na(promeasyscare), 
                                             promeasyscare + 1, NA)),
           promexistential = ifelse(!is.na(promexistential) & 
                                        promexistential == 15, NA,
                                      ifelse(!is.na(promexistential), 
                                             promexistential + 1, NA)),
           promnight = ifelse(!is.na(promnight) & 
                                  promnight == 15, NA,
                                      ifelse(!is.na(promnight), 
                                             promnight + 1, NA)),
           promnonsad = ifelse(!is.na(promnonsad) & 
                                   promnonsad == 15, NA,
                                      ifelse(!is.na(promnonsad), 
                                             promnonsad + 1, NA)), 
           promalone = ifelse(!is.na(promalone) & 
                                  promalone == 15, NA,
                                 ifelse(!is.na(promalone), 
                                        promalone + 1, NA)),
           promwrong = ifelse(!is.na(promwrong) & 
                                  promwrong == 15, NA,
                                 ifelse(!is.na(promwrong), 
                                        promwrong + 1, NA)),
           promright = ifelse(!is.na(promright) & 
                                  promright == 15, NA,
                                 ifelse(!is.na(promright), 
                                        promright + 1, NA)),
           promlonely = ifelse(!is.na(promlonely) & 
                                   promlonely == 15, NA,
                                 ifelse(!is.na(promlonely), 
                                        promlonely + 1, NA)),
           promsad = ifelse(!is.na(promsad) & 
                                promsad == 15, NA,
                                 ifelse(!is.na(promsad), 
                                        promsad + 1, NA)),
           promunhappy = ifelse(!is.na(promunhappy) & 
                                    promunhappy == 15, NA,
                                 ifelse(!is.na(promunhappy), 
                                        promunhappy + 1, NA)),
           promfun = ifelse(!is.na(promfun) & 
                                promfun == 15, NA,
                                 ifelse(!is.na(promfun), 
                                        promfun + 1, NA))) ->
    PROMISdata

# Make new columns with row means
PROMISdata$ped_anx_means <- rowMeans(subset(PROMISdata, 
                                            select=c("promawful", 
                                                     "promnervous", 
                                                     "promscared", 
                                                     "promworried", 
                                                     "promworryhome", 
                                                     "promeasyscare", 
                                                     "promexistential", 
                                                     "promnight")),
                                     na.rm=TRUE)
PROMISdata$ped_dep_means <- rowMeans(subset(PROMISdata, 
                                               select=c("promnonsad", 
                                                        "promalone", 
                                                        "promwrong", 
                                                        "promright", 
                                                        "promlonely", 
                                                        "promsad", 
                                                        "promunhappy", 
                                                        "promfun")),
                                     na.rm=TRUE)
PROMISdata$adult_anx_means <- rowMeans(subset(PROMISdata,
                                                 select=c("promisfear",
                                                          "promisanxious",
                                                          "promisworry",
                                                          "promisfocus",
                                                          "promisnervous",
                                                          "promisuneasy",
                                                          "promistense")),
                                       na.rm=TRUE)
PROMISdata$adult_dep_means <- rowMeans(subset(PROMISdata,
                                                 select=c("promisworthless",
                                                          "promisnothing",
                                                          "promissad",
                                                          "promisfailure",
                                                          "promisdepressed",
                                                          "promisunhappy",
                                                          "promishopeless",
                                                          "subject_promishelpless")),
                                       na.rm=TRUE)

# Make new columns with raw scores
PROMISdata$PROMIS_ped_anx <- PROMISdata$ped_anx_means*8
PROMISdata$PROMIS_ped_dep <- PROMISdata$ped_dep_means*8
PROMISdata$PROMIS_adult_anx <- PROMISdata$adult_anx_means*7
PROMISdata$PROMIS_adult_dep <- PROMISdata$adult_dep_means*8 

# Load conversion table
PROMISConversion <- read.csv("PROMISConversion.csv")

# Get conversion scores for each participant
PROMISdata %>%
    select(Subject.ID, Visit,
           PROMIS_ped_anx, PROMIS_ped_dep, PROMIS_adult_anx, 
           PROMIS_adult_dep) %>%
    pivot_longer(cols = c(PROMIS_ped_anx, PROMIS_ped_dep, PROMIS_adult_anx, 
                           PROMIS_adult_dep), 
                  names_to = "Scale", values_to = "Raw") %>%
    filter(!is.na(Raw)) %>%
    mutate(Raw = round(Raw)) %>%
    left_join(., PROMISConversion, by = c("Scale", "Raw")) %>%
    mutate(Scale = paste0(Scale, "_tscore")) %>%
    select(-Raw) %>%
    pivot_wider(names_from = Scale, values_from = Tscore) %>%
    full_join(PROMISdata, by = c("Subject.ID", "Visit")) %>% 
    arrange("Subject.ID", "Visit")->
    PROMISdata

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
write.csv(PROMISdata, file="PROMISdata.csv")

#### T-Tests ####

## All transgender youth on GAH
# Create dataframe with sample id and baseline PROMIS Anxiety values
baseline_PROMISanx <- PROMISdata[PROMISdata$Visit==1,
                                 c("Subject.ID",
                                   "PROMIS_ped_anx_tscore",
                                   "PROMIS_adult_anx_tscore")]
baseline_PROMISanx %>% unite(col="baseline_PROMISanx",PROMIS_ped_anx_tscore,
                             PROMIS_adult_anx_tscore,na.rm=T) ->
    baseline_PROMISanx
# Create dataframe with sample id and Visit 2 PROMIS Anxiety values
Visit2_PROMISanx <- PROMISdata[PROMISdata$Visit==2,
                               c("Subject.ID",
                                 "PROMIS_ped_anx_tscore",
                                 "PROMIS_adult_anx_tscore")]
Visit2_PROMISanx %>% unite(col="Visit2_PROMISanx",PROMIS_ped_anx_tscore,
                             PROMIS_adult_anx_tscore,na.rm=T) ->Visit2_PROMISanx
# Merge into one dataframe with 3 columns
PROMISanx_1_2 <- merge(baseline_PROMISanx, Visit2_PROMISanx, by.x="Subject.ID",
                       by.y="Subject.ID",all.x=F)
# Convert character columns to numeric
PROMISanx_1_2$baseline_PROMISanx <- 
    as.numeric(as.character(PROMISanx_1_2$baseline_PROMISanx))
PROMISanx_1_2$Visit2_PROMISanx <- 
    as.numeric(as.character(PROMISanx_1_2$Visit2_PROMISanx))
# Delete rows by name (delete those only on blockers)
PROMISanx_1_2 <- PROMISanx_1_2[!(row.names(PROMISanx_1_2) %in% 
                                     c("2","11","12")),]
saveRDS(PROMISanx_1_2,file="PROMISanx_1_2_GAH.RDS")
write.csv(PROMISanx_1_2,file="PROMISanx_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(PROMISanx_1_2$baseline_PROMISanx,
                 PROMISanx_1_2$Visit2_PROMISanx)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(PROMISanx_1_2$baseline_PROMISanx,
       PROMISanx_1_2$Visit2_PROMISanx,paired=T)

# Create dataframe with sample id and baseline PROMIS Depression values
baseline_PROMISdep <- PROMISdata[PROMISdata$Visit==1,
                                 c("Subject.ID",
                                   "PROMIS_ped_dep_tscore",
                                   "PROMIS_adult_dep_tscore")]
baseline_PROMISdep %>% unite(col="baseline_PROMISdep",PROMIS_ped_dep_tscore,
                             PROMIS_adult_dep_tscore,na.rm=T) ->
    baseline_PROMISdep
# Create dataframe with sample id and Visit 2 PROMIS Depression values
Visit2_PROMISdep <- PROMISdata[PROMISdata$Visit==2,c("Subject.ID",
                                                     "PROMIS_ped_dep_tscore",
                                                     "PROMIS_adult_dep_tscore")]
Visit2_PROMISdep %>% unite(col="Visit2_PROMISdep",PROMIS_ped_dep_tscore,
                           PROMIS_adult_dep_tscore,na.rm=T) ->Visit2_PROMISdep
# Merge into one dataframe with 3 columns
PROMISdep_1_2 <- merge(baseline_PROMISdep, Visit2_PROMISdep, by.x="Subject.ID",
                       by.y="Subject.ID",all.x=F)
# Convert character columns to numeric
PROMISdep_1_2$baseline_PROMISdep <- 
    as.numeric(as.character(PROMISdep_1_2$baseline_PROMISdep))
PROMISdep_1_2$Visit2_PROMISdep <- 
    as.numeric(as.character(PROMISdep_1_2$Visit2_PROMISdep))
# Delete rows by name (delete those only on blockers)
PROMISdep_1_2 <- PROMISdep_1_2[!(row.names(PROMISdep_1_2) %in% 
                                     c("2","11","12")),]
saveRDS(PROMISdep_1_2,file="PROMISdep_1_2_GAH.RDS")
write.csv(PROMISdep_1_2,file="PROMISdep_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(PROMISdep_1_2$baseline_PROMISdep,
                 PROMISdep_1_2$Visit2_PROMISdep)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(PROMISdep_1_2$baseline_PROMISdep,
       PROMISdep_1_2$Visit2_PROMISdep,paired=T)


## Transmasculine youth on GAH
# Create dataframe with sample id and baseline PROMIS Anxiety values
baseline_PROMISanx <- PROMISdata[PROMISdata$Visit==1,
                                 c("Subject.ID",
                                   "PROMIS_ped_anx_tscore",
                                   "PROMIS_adult_anx_tscore")]
baseline_PROMISanx %>% unite(col="baseline_PROMISanx",PROMIS_ped_anx_tscore,
                             PROMIS_adult_anx_tscore,na.rm=T)->
    baseline_PROMISanx
# Create dataframe with sample id and Visit 2 PROMIS Anxiety values
Visit2_PROMISanx <- PROMISdata[PROMISdata$Visit==2,c("Subject.ID",
                                                     "PROMIS_ped_anx_tscore",
                                                     "PROMIS_adult_anx_tscore")]
Visit2_PROMISanx %>% unite(col="Visit2_PROMISanx",PROMIS_ped_anx_tscore,
                           PROMIS_adult_anx_tscore,na.rm=T) ->Visit2_PROMISanx
# Merge into one dataframe with 3 columns
PROMISanx_1_2 <- merge(baseline_PROMISanx, Visit2_PROMISanx, by.x="Subject.ID",
                       by.y="Subject.ID",all.x=F)
# Convert character columns to numeric
PROMISanx_1_2$baseline_PROMISanx <- 
    as.numeric(as.character(PROMISanx_1_2$baseline_PROMISanx))
PROMISanx_1_2$Visit2_PROMISanx <- 
    as.numeric(as.character(PROMISanx_1_2$Visit2_PROMISanx))
# Delete rows by name (delete those only on blockers or assigned male at birth)
PROMISanx_1_2 <- PROMISanx_1_2[!(row.names(PROMISanx_1_2) %in% 
                                     c("1","2","9","11","12","15"
)),]
saveRDS(PROMISanx_1_2,file="PROMISanx_1_2_GAH.RDS")
write.csv(PROMISanx_1_2,file="PROMISanx_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(PROMISanx_1_2$baseline_PROMISanx,
                 PROMISanx_1_2$Visit2_PROMISanx)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(PROMISanx_1_2$baseline_PROMISanx,
       PROMISanx_1_2$Visit2_PROMISanx,paired=T)

# Create dataframe with sample id and baseline PROMIS Depression values
baseline_PROMISdep <- PROMISdata[PROMISdata$Visit==1,
                                 c("Subject.ID",
                                   "PROMIS_ped_dep_tscore",
                                   "PROMIS_adult_dep_tscore")]
baseline_PROMISdep %>% unite(col="baseline_PROMISdep",PROMIS_ped_dep_tscore,
                             PROMIS_adult_dep_tscore,
                             na.rm=T)-> baseline_PROMISdep
# Create dataframe with sample id and Visit 2 PROMIS Depression values
Visit2_PROMISdep <- PROMISdata[PROMISdata$Visit==2,c("Subject.ID",
                                                     "PROMIS_ped_dep_tscore",
                                                     "PROMIS_adult_dep_tscore")]
Visit2_PROMISdep %>% unite(col="Visit2_PROMISdep",PROMIS_ped_dep_tscore,
                           PROMIS_adult_dep_tscore,na.rm=T) ->Visit2_PROMISdep
# Merge into one dataframe with 3 columns
PROMISdep_1_2 <- merge(baseline_PROMISdep, Visit2_PROMISdep, by.x="Subject.ID",
                       by.y="Subject.ID",all.x=F)
# Convert character columns to numeric
PROMISdep_1_2$baseline_PROMISdep <- 
    as.numeric(as.character(PROMISdep_1_2$baseline_PROMISdep))
PROMISdep_1_2$Visit2_PROMISdep <- 
    as.numeric(as.character(PROMISdep_1_2$Visit2_PROMISdep))
# Delete rows by name (delete those only on blockers or assigned male at birth)
PROMISdep_1_2 <- PROMISdep_1_2[!(row.names(PROMISdep_1_2) %in% 
                                     c("1","2","9","11","12","15"
)),]
saveRDS(PROMISdep_1_2,file="PROMISdep_1_2_GAH.RDS")
write.csv(PROMISdep_1_2,file="PROMISdep_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(PROMISdep_1_2$baseline_PROMISdep,
                 PROMISdep_1_2$Visit2_PROMISdep)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(PROMISdep_1_2$baseline_PROMISdep,
       PROMISdep_1_2$Visit2_PROMISdep,paired=T)