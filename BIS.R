# Body Image Scale

# Load packages
library(dplyr)
library(tidyr)
library(piggyback)
library(readr)
library(base)
library(stringr)

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data")

# Load data
BISdata <- read.csv("BISdata.csv")

# Create object for desire to change body part
BISdata %>% select(Subject.ID:Visit, subject_ynnose:subject_ynstature) %>% 
    replace(., .==14 | .==15,NA) ->
    BISdata_change

# Make new columns with row means for each category of body part
BISdata_change$primary <- rowMeans(subset(BISdata_change, 
                                        select=c("subject_ynpenis", 
                                                 "subject_ynscrotum", 
                                                 "subject_yntesticles", 
                                                 "subject_ynfacialhair", 
                                                 "subject_ynbodyhair", 
                                                 "subject_ynbreasts", 
                                                 "subject_ynvagina", 
                                                 "subject_ynclitoris",
                                                 "subject_ynovaries",
                                                 "subject_ynvoice",
                                                 "subject_ynchest")),
                                 na.rm=TRUE)
BISdata_change$secondary <- rowMeans(subset(BISdata_change, 
                                          select=c("subject_ynhips", 
                                                   "subject_ynfigure", 
                                                   "subject_ynwaist", 
                                                   "subject_ynarms", 
                                                   "subject_ynbutt", 
                                                   "subject_ynbiceps", 
                                                   "subject_ynappearance", 
                                                   "subject_ynstature",
                                                   "subject_ynmuscles",
                                                   "subject_ynweight",
                                                   "subject_ynthighs",
                                                   "subject_ynhair")),
                                   na.rm=TRUE)
BISdata_change$neutral <- rowMeans(subset(BISdata_change, 
                                        select=c("subject_ynnose", 
                                                 "subject_ynshoulders", 
                                                 "subject_ynchin", 
                                                 "subject_yncalves", 
                                                 "subject_ynhands", 
                                                 "subject_ynadamsapple", 
                                                 "subject_yneyebrows", 
                                                 "subject_ynface",
                                                 "subject_ynfeet",
                                                 "subject_ynheight")),
                                 na.rm=TRUE)



# Create object for participants' dissatisfaction with body parts
BISdata %>%
    select(Subject.ID, Visit, contains("___")) %>% 
    # Eliminate columns indicating skip columns saying you don't have body part
    select(-contains("___15"), -contains("___0")) %>% 
    gather(Variable, Value, -Subject.ID, -Visit) %>%
    # Keep only instances that are CHOSEN
    # Eliminate all non-1
    filter(!is.na(Value)) %>%
    filter(Value != 0) %>%
    select(-Value) %>%
    # Extract needed information from the VARIABLE NAME
    # Eliminate the start of variable name
    mutate(Variable = str_remove(Variable, "subject_")) %>% 
    # Separate out body part from response
    separate(Variable, into = c("BodyPart", "OptionChosen"), sep = "___") %>% 
    # Turn chosen response into a number
    mutate(OptionChosen = as.numeric(OptionChosen)) %>%
    # In some cases, participants chose multiple options 
    # (e.g., 10001 on visit 2 chose 2 and 3 about their butt)
    # In these instances, we take the mean
    group_by(Subject.ID, Visit, BodyPart) %>% 
    summarise(OptionChosen = mean(OptionChosen),
              .groups = "drop") %>%
    pivot_wider(names_from = BodyPart, values_from = OptionChosen)->
    BISdata_feel

# Make new columns with row means for each category of body part
BISdata_feel$primary <- rowMeans(subset(BISdata_feel, 
                                        select=c("penis", 
                                                 "scrotum", 
                                                 "testicles", 
                                                 "facialhair", 
                                                 "bodyhair", 
                                                 "breasts", 
                                                 "vagina", 
                                                 "clitoris",
                                                 "ovaries",
                                                 "voice",
                                                 "chest")),
                                 na.rm=TRUE)
BISdata_feel$secondary <- rowMeans(subset(BISdata_feel, 
                                        select=c("hips", 
                                                 "figure", 
                                                 "waist", 
                                                 "arms", 
                                                 "butt", 
                                                 "biceps", 
                                                 "appearance", 
                                                 "stature",
                                                 "muscles",
                                                 "weight",
                                                 "thighs",
                                                 "hair")),
                                 na.rm=TRUE)
BISdata_feel$neutral <- rowMeans(subset(BISdata_feel, 
                                          select=c("nose", 
                                                   "shoulders", 
                                                   "chin", 
                                                   "calves", 
                                                   "hands", 
                                                   "adamsapple", 
                                                   "eyebrows", 
                                                   "face",
                                                   "feet",
                                                   "height")),
                                   na.rm=TRUE)

# Set working directory
setwd("/Users/gillianhilscher/Downloads/Princeton_Spring_2023/Thesis/DATA/Psych Data/03-24-23")
write.csv(BISdata_feel,file="BISdata_feel.csv")
write.csv(BISdata_change,file="BISdata_change.csv")

#### T-Tests ####

## All transgender youth on GAH
# Create dataframe with sample id and baseline BIS feel primary values
baseline_BISfeel_primary <- BISdata_feel[BISdata_feel$Visit==1,c("Subject.ID","primary")]
names(baseline_BISfeel_primary)[2]<-"baseline_BISfeel_primary"
# Create dataframe with sample id and Visit 2 BIS feel secondary values
Visit2_BISfeel_primary <- BISdata_feel[BISdata_feel$Visit==2,c("Subject.ID","primary")]
names(Visit2_BISfeel_primary)[2]<-"Visit2_BISfeel_primary"
# Merge into one dataframe with 3 columns
BISfeel_primary_1_2 <- merge(baseline_BISfeel_primary, 
                             Visit2_BISfeel_primary, by.x="Subject.ID",
                             by.y="Subject.ID",all.x=F)
# Delete rows by name (delete 10022 because answered NA to all primary questions)
BISfeel_primary_1_2 <- BISfeel_primary_1_2[!(row.names(BISfeel_primary_1_2) %in%
                                                 c("15")),]
# Delete rows by name (delete those only on blockers)
BISfeel_primary_1_2 <- BISfeel_primary_1_2[!(row.names(BISfeel_primary_1_2)%in%
                                                 c("2","11","12")),]
saveRDS(BISfeel_primary_1_2,file="BISfeel_primary_1_2_GAH.RDS")
write.csv(BISfeel_primary_1_2,file="BISfeel_primary_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISfeel_primary_1_2$baseline_BISfeel_primary,
                 BISfeel_primary_1_2$Visit2_BISfeel_primary)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISfeel_primary_1_2$baseline_BISfeel_primary,
       BISfeel_primary_1_2$Visit2_BISfeel_primary,paired=T)

# Create dataframe with sample id and baseline BIS feel secondary values
baseline_BISfeel_secondary <- BISdata_feel[BISdata_feel$Visit==1,c("Subject.ID","secondary")]
names(baseline_BISfeel_secondary)[2]<-"baseline_BISfeel_secondary"
# Create dataframe with sample id and Visit 2 BIS feel secondary values
Visit2_BISfeel_secondary <- BISdata_feel[BISdata_feel$Visit==2,c("Subject.ID","secondary")]
names(Visit2_BISfeel_secondary)[2]<-"Visit2_BISfeel_secondary"
# Merge into one dataframe with 3 columns
BISfeel_secondary_1_2 <- merge(baseline_BISfeel_secondary, 
                               Visit2_BISfeel_secondary, by.x="Subject.ID",
                               by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
BISfeel_secondary_1_2<-BISfeel_secondary_1_2[!(row.names(BISfeel_secondary_1_2)%in%
                                                   c("2","11","12")),]
saveRDS(BISfeel_secondary_1_2,file="BISfeel_secondary_1_2_GAH.RDS")
write.csv(BISfeel_secondary_1_2,file="BISfeel_secondary_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISfeel_secondary_1_2$baseline_BISfeel_secondary,
                 BISfeel_secondary_1_2$Visit2_BISfeel_secondary)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISfeel_secondary_1_2$baseline_BISfeel_secondary,
       BISfeel_secondary_1_2$Visit2_BISfeel_secondary,paired=T)

# Create dataframe with sample id and baseline BIS feel neutral values
baseline_BISfeel_neutral <- BISdata_feel[BISdata_feel$Visit==1,c("Subject.ID","neutral")]
names(baseline_BISfeel_neutral)[2]<-"baseline_BISfeel_neutral"
# Create dataframe with sample id and Visit 2 BIS feel neutral values
Visit2_BISfeel_neutral <- BISdata_feel[BISdata_feel$Visit==2,c("Subject.ID","neutral")]
names(Visit2_BISfeel_neutral)[2]<-"Visit2_BISfeel_neutral"
# Merge into one dataframe with 3 columns
BISfeel_neutral_1_2 <- merge(baseline_BISfeel_neutral, Visit2_BISfeel_neutral, 
                             by.x="Subject.ID",by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
BISfeel_neutral_1_2 <- BISfeel_neutral_1_2[!(row.names(BISfeel_neutral_1_2)%in%
                                                 c("2","11","12")),]
saveRDS(BISfeel_neutral_1_2,file="BISfeel_neutral_1_2_GAH.RDS")
write.csv(BISfeel_neutral_1_2,file="BISfeel_neutral_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISfeel_neutral_1_2$baseline_BISfeel_neutral,
                 BISfeel_neutral_1_2$Visit2_BISfeel_neutral)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISfeel_neutral_1_2$baseline_BISfeel_neutral,
       BISfeel_neutral_1_2$Visit2_BISfeel_neutral,paired=T)

# Create dataframe with sample id and baseline BIS change primary values
baseline_BISchange_primary <- BISdata_change[BISdata_change$Visit==1,c("Subject.ID","primary")]
names(baseline_BISchange_primary)[2]<-"baseline_BISchange_primary"
# Create dataframe with sample id and Visit 2 BIS change secondary values
Visit2_BISchange_primary <- BISdata_change[BISdata_change$Visit==2,c("Subject.ID","primary")]
names(Visit2_BISchange_primary)[2]<-"Visit2_BISchange_primary"
# Merge into one dataframe with 3 columns
BISchange_primary_1_2 <- merge(baseline_BISchange_primary, 
                               Visit2_BISchange_primary, by.x="Subject.ID",
                               by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
BISchange_primary_1_2<-BISchange_primary_1_2[!(row.names(BISchange_primary_1_2)%in%
                                                   c("2","11","12")),]
saveRDS(BISchange_primary_1_2,file="BISchange_primary_1_2_GAH.RDS")
write.csv(BISchange_primary_1_2,file="BISchange_primary_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISchange_primary_1_2$baseline_BISchange_primary,
                 BISchange_primary_1_2$Visit2_BISchange_primary)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISchange_primary_1_2$baseline_BISchange_primary,
       BISchange_primary_1_2$Visit2_BISchange_primary,paired=T)

# Create dataframe with sample id and baseline BIS change secondary values
baseline_BISchange_secondary <- BISdata_change[BISdata_change$Visit==1,c("Subject.ID","secondary")]
names(baseline_BISchange_secondary)[2]<-"baseline_BISchange_secondary"
# Create dataframe with sample id and Visit 2 BIS change secondary values
Visit2_BISchange_secondary <- BISdata_change[BISdata_change$Visit==2,c("Subject.ID","secondary")]
names(Visit2_BISchange_secondary)[2]<-"Visit2_BISchange_secondary"
# Merge into one dataframe with 3 columns
BISchange_secondary_1_2 <- merge(baseline_BISchange_secondary, 
                                 Visit2_BISchange_secondary, by.x="Subject.ID",
                                 by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
BISchange_secondary_1_2<-BISchange_secondary_1_2[!(row.names(BISchange_secondary_1_2)%in%
                                                       c("2","11","12")),]
saveRDS(BISchange_secondary_1_2,file="BISchange_secondary_1_2_GAH.RDS")
write.csv(BISchange_secondary_1_2,file="BISchange_secondary_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISchange_secondary_1_2$baseline_BISchange_secondary,
                 BISchange_secondary_1_2$Visit2_BISchange_secondary)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISchange_secondary_1_2$baseline_BISchange_secondary,
       BISchange_secondary_1_2$Visit2_BISchange_secondary,paired=T)

# Create dataframe with sample id and baseline BIS change neutral values
baseline_BISchange_neutral <- BISdata_change[BISdata_change$Visit==1,c("Subject.ID","neutral")]
names(baseline_BISchange_neutral)[2]<-"baseline_BISchange_neutral"
# Create dataframe with sample id and Visit 2 BIS change neutral values
Visit2_BISchange_neutral <- BISdata_change[BISdata_change$Visit==2,c("Subject.ID","neutral")]
names(Visit2_BISchange_neutral)[2]<-"Visit2_BISchange_neutral"
# Merge into one dataframe with 3 columns
BISchange_neutral_1_2 <- merge(baseline_BISchange_neutral, 
                               Visit2_BISchange_neutral, by.x="Subject.ID",
                               by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers)
BISchange_neutral_1_2<-BISchange_neutral_1_2[!(row.names(BISchange_neutral_1_2)%in%
                                                   c("2","11","12")),]
saveRDS(BISchange_neutral_1_2,file="BISchange_neutral_1_2_GAH.RDS")
write.csv(BISchange_neutral_1_2,file="BISchange_neutral_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISchange_neutral_1_2$baseline_BISchange_neutral,
                 BISchange_neutral_1_2$Visit2_BISchange_neutral)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISchange_neutral_1_2$baseline_BISchange_neutral,
       BISchange_neutral_1_2$Visit2_BISchange_neutral,paired=T)


## Transmasculine youth on GAH
# Create dataframe with sample id and baseline BIS feel primary values
baseline_BISfeel_primary <- BISdata_feel[BISdata_feel$Visit==1,c("Subject.ID","primary")]
names(baseline_BISfeel_primary)[2]<-"baseline_BISfeel_primary"
# Create dataframe with sample id and Visit 2 BIS feel secondary values
Visit2_BISfeel_primary <- BISdata_feel[BISdata_feel$Visit==2,c("Subject.ID","primary")]
names(Visit2_BISfeel_primary)[2]<-"Visit2_BISfeel_primary"
# Merge into one dataframe with 3 columns
BISfeel_primary_1_2 <- merge(baseline_BISfeel_primary, Visit2_BISfeel_primary, 
                             by.x="Subject.ID",by.y="Subject.ID",all.x=F)
# Delete rows by name (delete 10022 because answered NA to all primary questions)
BISfeel_primary_1_2 <- BISfeel_primary_1_2[!(row.names(BISfeel_primary_1_2)%in%
                                                 c("15")),]
# Delete rows by name (delete those only on blockers or assigned male at birth)
BISfeel_primary_1_2 <- BISfeel_primary_1_2[!(row.names(BISfeel_primary_1_2)%in%
                                                 c("1","2","9","11","12","16")),]
saveRDS(BISfeel_primary_1_2,file="BISfeel_primary_1_2_GAH.RDS")
write.csv(BISfeel_primary_1_2,file="BISfeel_primary_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISfeel_primary_1_2$baseline_BISfeel_primary,
                 BISfeel_primary_1_2$Visit2_BISfeel_primary)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISfeel_primary_1_2$baseline_BISfeel_primary,
       BISfeel_primary_1_2$Visit2_BISfeel_primary,paired=T)


# Create dataframe with sample id and baseline BIS feel secondary values
baseline_BISfeel_secondary <- BISdata_feel[BISdata_feel$Visit==1,
                                           c("Subject.ID","secondary")]
names(baseline_BISfeel_secondary)[2]<-"baseline_BISfeel_secondary"
# Create dataframe with sample id and Visit 2 BIS feel secondary values
Visit2_BISfeel_secondary <- BISdata_feel[BISdata_feel$Visit==2,
                                         c("Subject.ID","secondary")]
names(Visit2_BISfeel_secondary)[2]<-"Visit2_BISfeel_secondary"
# Merge into one dataframe with 3 columns
BISfeel_secondary_1_2 <- merge(baseline_BISfeel_secondary, 
                               Visit2_BISfeel_secondary, 
                               by.x="Subject.ID",by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
BISfeel_secondary_1_2<-BISfeel_secondary_1_2[!(row.names(BISfeel_secondary_1_2)%in%
                                                   c("1","2","9","11","12","16")),]
saveRDS(BISfeel_secondary_1_2,file="BISfeel_secondary_1_2_GAH.RDS")
write.csv(BISfeel_secondary_1_2,file="BISfeel_secondary_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISfeel_secondary_1_2$baseline_BISfeel_secondary,
                 BISfeel_secondary_1_2$Visit2_BISfeel_secondary)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISfeel_secondary_1_2$baseline_BISfeel_secondary,
       BISfeel_secondary_1_2$Visit2_BISfeel_secondary,paired=T)

# Create dataframe with sample id and baseline BIS feel neutral values
baseline_BISfeel_neutral <- BISdata_feel[BISdata_feel$Visit==1,
                                         c("Subject.ID","neutral")]
names(baseline_BISfeel_neutral)[2]<-"baseline_BISfeel_neutral"
# Create dataframe with sample id and Visit 2 BIS feel neutral values
Visit2_BISfeel_neutral <- BISdata_feel[BISdata_feel$Visit==2,
                                       c("Subject.ID","neutral")]
names(Visit2_BISfeel_neutral)[2]<-"Visit2_BISfeel_neutral"
# Merge into one dataframe with 3 columns
BISfeel_neutral_1_2 <- merge(baseline_BISfeel_neutral, Visit2_BISfeel_neutral, 
                             by.x="Subject.ID",by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
BISfeel_neutral_1_2 <- BISfeel_neutral_1_2[!(row.names(BISfeel_neutral_1_2)%in%
                                                 c("1","2","9","11","12","16")),]
saveRDS(BISfeel_neutral_1_2,file="BISfeel_neutral_1_2_GAH.RDS")
write.csv(BISfeel_neutral_1_2,file="BISfeel_neutral_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISfeel_neutral_1_2$baseline_BISfeel_neutral,
                 BISfeel_neutral_1_2$Visit2_BISfeel_neutral)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISfeel_neutral_1_2$baseline_BISfeel_neutral,
       BISfeel_neutral_1_2$Visit2_BISfeel_neutral,paired=T)

# Create dataframe with sample id and baseline BIS change primary values
baseline_BISchange_primary <- BISdata_change[BISdata_change$Visit==1,
                                             c("Subject.ID","primary")]
names(baseline_BISchange_primary)[2]<-"baseline_BISchange_primary"
# Create dataframe with sample id and Visit 2 BIS change secondary values
Visit2_BISchange_primary <- BISdata_change[BISdata_change$Visit==2,
                                           c("Subject.ID","primary")]
names(Visit2_BISchange_primary)[2]<-"Visit2_BISchange_primary"
# Merge into one dataframe with 3 columns
BISchange_primary_1_2 <- merge(baseline_BISchange_primary, 
                               Visit2_BISchange_primary, by.x="Subject.ID",
                               by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
BISchange_primary_1_2<-BISchange_primary_1_2[!(row.names(BISchange_primary_1_2)%in%
                                                     c("1","2","9","11","12","16")),]
saveRDS(BISchange_primary_1_2,file="BISchange_primary_1_2_GAH.RDS")
write.csv(BISchange_primary_1_2,file="BISchange_primary_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISchange_primary_1_2$baseline_BISchange_primary,
                 BISchange_primary_1_2$Visit2_BISchange_primary)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISchange_primary_1_2$baseline_BISchange_primary,
       BISchange_primary_1_2$Visit2_BISchange_primary,paired=T)

# Create dataframe with sample id and baseline BIS change secondary values
baseline_BISchange_secondary <- BISdata_change[BISdata_change$Visit==1,c("Subject.ID","secondary")]
names(baseline_BISchange_secondary)[2]<-"baseline_BISchange_secondary"
# Create dataframe with sample id and Visit 2 BIS change secondary values
Visit2_BISchange_secondary <- BISdata_change[BISdata_change$Visit==2,c("Subject.ID","secondary")]
names(Visit2_BISchange_secondary)[2]<-"Visit2_BISchange_secondary"
# Merge into one dataframe with 3 columns
BISchange_secondary_1_2 <- merge(baseline_BISchange_secondary, 
                                 Visit2_BISchange_secondary, 
                                 by.x="Subject.ID",by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
BISchange_secondary_1_2<-BISchange_secondary_1_2[!(row.names(BISchange_secondary_1_2)%in%
                                                       c("1","2","9","11","12","16")),]
saveRDS(BISchange_secondary_1_2,file="BISchange_secondary_1_2_GAH.RDS")
write.csv(BISchange_secondary_1_2,file="BISchange_secondary_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISchange_secondary_1_2$baseline_BISchange_secondary,
                 BISchange_secondary_1_2$Visit2_BISchange_secondary)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISchange_secondary_1_2$baseline_BISchange_secondary,
       BISchange_secondary_1_2$Visit2_BISchange_secondary,paired=T)

# Create dataframe with sample id and baseline BIS change neutral values
baseline_BISchange_neutral <- BISdata_change[BISdata_change$Visit==1,c("Subject.ID","neutral")]
names(baseline_BISchange_neutral)[2]<-"baseline_BISchange_neutral"
# Create dataframe with sample id and Visit 2 BIS change neutral values
Visit2_BISchange_neutral <- BISdata_change[BISdata_change$Visit==2,c("Subject.ID","neutral")]
names(Visit2_BISchange_neutral)[2]<-"Visit2_BISchange_neutral"
# Merge into one dataframe with 3 columns
BISchange_neutral_1_2<-merge(baseline_BISchange_neutral, 
                             Visit2_BISchange_neutral, by.x="Subject.ID",
                             by.y="Subject.ID",all.x=F)
# Delete rows by name (delete those only on blockers or assigned male at birth)
BISchange_neutral_1_2<-BISchange_neutral_1_2[!(row.names(BISchange_neutral_1_2)%in%
                                                     c("1","2","9","11","12","16")),]
saveRDS(BISchange_neutral_1_2,file="BISchange_neutral_1_2_GAH.RDS")
write.csv(BISchange_neutral_1_2,file="BISchange_neutral_1_2_GAH.csv")
# Perform Shapiro test to see distribution
shapiro.test(x=c(BISchange_neutral_1_2$baseline_BISchange_neutral,
                 BISchange_neutral_1_2$Visit2_BISchange_neutral)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(BISchange_neutral_1_2$baseline_BISchange_neutral,
       BISchange_neutral_1_2$Visit2_BISchange_neutral,paired=T)