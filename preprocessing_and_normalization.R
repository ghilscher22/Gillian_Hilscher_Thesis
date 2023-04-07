# DNAm Pre-Processing and Normalization

# Load packages
library(minfi)
library(ENmix)
library(dplyr)
library(missMethyl)

#####  Normalization/Quality Control of Trans Youth Data in ENmix ####
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/All idats/")

# Make dataframe of IDAT files
df <- data.frame(basename = list.files(pattern="idat"))
head(df)
# Remove Grn/Red labels
df$basename <- gsub("_(Grn|Red).idat","",df$basename)
# Remove duplicates
df <- df[!duplicated(df$basename),]
# Create dataframe
df <- as.data.frame(df)
# Change column name to "basename"
names(df) <- "basename"

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

# Create RGChannelSetExtended
rg<-read.metharray(basenames=df$basename,extended=T)
saveRDS(rg,file="G954_EPIC_Data_0217rgsetext.RDS")

# Attach metadata for later use
baseDir <- file.path("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-20-23")
G954sheet <- read.metharray.sheet(base=baseDir,
                                  pattern = "G954_Sample_Sheet.csv")
rownames(G954sheet)=G954sheet$Sample_Name
colData(rg)=as(G954sheet[colnames(rg),],"DataFrame")
G954sheet %>% as.data.frame() %>% 
    arrange(Sample_Name)->G954sheet


# Signal intensity quality control
qc<-QCinfo(rg)
str(qc)


# Estimate cell type proportions using ENmix
celltype_ENmix<-estimateCellProp(userdata=rg,
                                 refdata="FlowSorted.Blood.EPIC",
                                 nonnegative = TRUE,normalize=TRUE)
# Delete rows by name (delete sample removed by qc)
celltype_ENmix <- celltype_ENmix[!(row.names(celltype_ENmix) %in% 
                                       c("206839580077_R08C01")),]
write.csv(celltype_ENmix,file="G954_EPIC_Data_0217celltype_ENmix.csv",
          row.names=FALSE)
saveRDS(celltype_ENmix,file="G954_EPIC_Data_0217celltype_ENmix.RDS")


# Perform background correction and dye bias correction using internal control 
# probes in ENmix
rg_betas <- preprocessENmix(rg, QCinfo = qc, nCores = 1) %>% 
    rcp(nbthre = 4) # probe-type bias adjustment
saveRDS(rg_betas,file="G954_EPIC_Data_0217rg_betas.RDS")
write.csv(rg_betas,file="G954_EPIC_Data_0217rg_betas.csv")



## Sex prediction in ENmix
# Create rgSet; plot and save RDS
path <- file.path("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/All idats/")
rgSet <- readidat(path = path,recursive = TRUE)
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23/PlotCtrl")
plotCtrl(rgSet)
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
saveRDS(rgSet,file="G954_EPIC_Data_0217rgSet.RDS")

# Predict sex in ENmix
rg_predSex <- predSex(rgSet)
write.csv(rg_predSex,file="G954_EPIC_Data_0217rg_predSex.csv")
saveRDS(rg_predSex,file="G954_EPIC_Data_0217rg_predSex.RDS")


## Methylation age calculation in ENmix
mage <- methyAge(rg_betas)
write.csv(mage,file="G954_EPIC_Data_0217mage.csv")
saveRDS(mage,file="G954_EPIC_Data_0217mage.RDS")


## Merge cell type, methylation age, and metadata for later use
# Merge methylation age and cell type data
mage_celltype <- merge(mage, celltype, by.x="SampleID",by.y="SampleID")
# Merge with metadata
G954_mage_ageAcc<-merge(mage_celltype,G954sheet,by.x="Sample_Name",
                        by.y="Sample_Name")
saveRDS(G954_mage_ageAcc,file="G954_mage_ageAcc.RDS")
write.csv(G954_mage_ageAcc,file="G954_mage_ageAcc.csv")


#####  Normalization/Quality Control of Trans Adult Data in ENmix ####

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/GSE173717_RAW/")

# Make dataframe of IDAT files
df <- data.frame(basename = list.files(pattern="idat"))
# Remove Grn/Red labels
df$basename <- gsub("_(Grn|Red).idat","",df$basename)
# Remove duplicates
df <- df[!duplicated(df$basename),]
# Create dataframe
df <- as.data.frame(df)
# Change column name to "basename"
names(df) <- "basename"

# Create RGChannelSetExtended
rg<-read.metharray(basenames=df$basename,extended=T,force=T)

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/03-28-23")
saveRDS(rg,file="GSE173717rgsetext.RDS")

# Attach metadata for later use
baseDir <- file.path("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/03-28-23")
sheet <- read.metharray.sheet(base=baseDir,pattern = "GSE173717_Sample_Sheet.csv")
rownames(sheet)=sheet$Sample_Name
colData(rg)=as(sheet[colnames(rg),],"DataFrame")
sheet %>% as.data.frame() %>% 
    arrange(Sample_Name)->sheet
saveRDS(sheet,file="GSE173717_Sample_Sheet.RDS")


# Signal intensity quality control
qc<-QCinfo(rg)
str(qc)


# Estimate cell type proportions with ENmix
celltype_ENmix<-estimateCellProp(userdata=rg,
                                 refdata="FlowSorted.Blood.EPIC",
                                 nonnegative = TRUE,normalize=TRUE)
write.csv(celltype_ENmix,file="GSE173717celltype_ENmix.csv",row.names=FALSE)
saveRDS(celltype_ENmix,file="GSE173717celltype_ENmix.RDS")

# Perform background correction and dye bias correction using internal control 
# probes in ENmix
rg_betas <- preprocessENmix(rg, QCinfo = qc, nCores = 1) %>% 
    rcp(nbthre = 4) # probe-type bias adjustment
saveRDS(rg_betas,file="GSE173717rg_betas.RDS")
write.csv(rg_betas,file="GSE173717rg_betas.csv")



## Sex prediction in ENmix
# Create rgSet; plot and save RDS
path <- file.path("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/GSE173717_RAW/")
rgSet <- readidat(path = path,recursive = TRUE)
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/03-28-23/PlotCtrl")
plotCtrl(rgSet)
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/03-28-23")
saveRDS(rgSet,file="GSE173717rgSet.RDS")
class(rgSet)

# Predict sex in ENmix
rg_predSex <- predSex(rgSet)
write.csv(rg_predSex,file="GSE173717rg_predSex.csv")
saveRDS(rg_predSex,file="GSE173717rg_predSex.RDS")


## Methylation age calculation in ENmix
mage <- methyAge(rg_betas)
write.csv(mage,file="GSE173717mage.csv")
saveRDS(mage,file="GSE173717mage.RDS")