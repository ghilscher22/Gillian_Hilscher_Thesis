# EWAS/Manhattan Plot (adapted from Jaeyoon Cha's thesis)

# Load packages
library(dplyr)
library(missMethyl)
library(ENmix)
library(limma)
library(ggplot2)
library(readr)
library(piggyback)
library(tidyverse)
library(qqman)
library(edgeR)

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

###### Differential Methylation in All Transgender Youth on GAH ######
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

# Create metadata object for Visit 1 and Visit 2: GAH only
G954_all<-G954_GAH[!G954_GAH$Visit %in% 3,]
G954_all<-G954_all %>% arrange(Sample_Name)
# On testosterone
G954_allTB<-G954_GAH[!G954_GAH$Visit %in% 3,]
G954_allTB<-G954_allTB[!G954_allTB$Female %in% 0,]
G954_allTB<-G954_allTB %>% arrange(Sample_Name)

# Create metadata object for Visit 1 and Visit 3: GAH only
G954_13<-G954_GAH %>% subset(.,Sample_Name %in%
                                   c("206842050024_R03C01", #10001
                                     "206842050024_R02C01", #10001_3
                                     "206839580077_R05C01", #10004
                                     "206839580077_R04C01", #10004_3
                                     "206839580077_R01C01", #10005
                                     "206839580077_R03C01")) #10005_3
G954_13<-G954_13 %>% arrange(Sample_Name)
# On testosterone
G954_13TB<-G954_GAH%>% subset(.,Sample_Name %in%
                                   c("206842050024_R03C01", #10001
                                     "206842050024_R02C01", #10001_3
                                     "206839580077_R05C01", #10004
                                     "206839580077_R04C01", #10004_3
                                     "206839580077_R01C01", #10005
                                     "206839580077_R03C01")) #10005_3
G954_13TB<-G954_13TB[!G954_13TB$Female %in% 0,]
G954_13TB<-G954_13TB %>% arrange(Sample_Name)

# Create vector of Sample_Names from Visit 1/2 samples for all GAH and 
# testosterone
string_use <- G954_all[,"Sample_Name"]
string_useTB <- G954_allTB[,"Sample_Name"]
# Create vector of Sample_Names for paired Visit 1/3 data
string_13 <- G954_13[,"Sample_Name"]
string_13TB<-G954_13TB[,"Sample_Name"]

# Load cell type proportion data
#celltype<-readRDS("G954_EPIC_Data_0217celltype_ENmix.RDS")
# Keep only cell types with Sample_Name in string_use
#celltype<-celltype[celltype$Sample_Name %in% string_use,] 
# Arrange celltype and metadata by Sample_Name
#celltype<-celltype %>% arrange(Sample_Name)

# Select M-value columns for all GAH visits 1 and 2
Mval_all12<-Mval[,(colnames(Mval) %in% string_use)] 
# On testosterone
Mval_all12TB<-Mval[,(colnames(Mval) %in% string_useTB)] 
# Select M-value columns for all GAH visits 1 and 3
Mval_all13<-Mval[,(colnames(Mval) %in% string_13)] 
# On testosterone
Mval_all13TB<-Mval[,(colnames(Mval) %in% string_13TB)] 

# Load manifest data (from Benjamin Jacobson's thesis)
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/zhou")
man <- readRDS("EPIC.hg19.manifest.rds")
man <- as(man, "data.frame")
man$Probe_ID <- rownames(man)
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

# Exclude bad probes and sex chromosomes (from Benjamin Jacobson's thesis)
manSub <- man[!man$MASK_general,]
manSub <- manSub[manSub$seqnames != "chrX",]
manSub <- manSub[manSub$seqnames != "chrY",]
manSub <- manSub[manSub$seqnames != "chrM",]

# Exclude all bad probes, SNPs, and sex chromosomes from M-value data (from 
# Benjamin Jacobson's thesis)
# Visit 1 and 2 data
Mval_all12_clean<-Mval_all12[rownames(Mval_all12) %in% manSub$Probe_ID,] 
# On testosterone
Mval_all12TB_clean<-Mval_all12TB[rownames(Mval_all12TB) %in% manSub$Probe_ID,]
# Visit 1 and Visit 3 data
Mval_all13_clean<-Mval_all13[rownames(Mval_all13) %in% manSub$Probe_ID,] 
# On testosterone
Mval_all13TB_clean<-Mval_all13TB[rownames(Mval_all13TB) %in% manSub$Probe_ID,]

# Remove probes that don't change: Visit 1 and 2
std<-apply(Mval_all12_clean,1,sd)
cutoff<-summary(std)[2]/2
keep<-std>cutoff
Mval_all122<-Mval_all12_clean[keep,]
# On testosterone
std<-apply(Mval_all12TB_clean,1,sd)
cutoff<-summary(std)[2]/2
keep<-std>cutoff
Mval_all122TB<-Mval_all12TB_clean[keep,]

# Remove probes that don't change: Visit 1 and 3
std<-apply(Mval_all13_clean,1,sd)
cutoff<-summary(std)[2]/2
keep<-std>cutoff
Mval_all133<-Mval_all13_clean[keep,]
# On testosterone
std<-apply(Mval_all13TB_clean,1,sd)
cutoff<-summary(std)[2]/2
keep<-std>cutoff
Mval_all133TB<-Mval_all13TB_clean[keep,]

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/EWAS")

############# Visit 1 and 2 #############

# Create column in metadata that combines Subject.ID with Visit
G954_all<-unite(G954_all, col="ID.Visit", c("Subject.ID","Visit"), sep="_",
                remove=F)
# Rename Mval_all122 columns with ID.Visit matches
colnames(Mval_all122)<-G954_all$ID.Visit[match(colnames(Mval_all122), 
                                               G954_all$Sample_Name)]

# Place Mval columns in numerical/Visit order
col_order <- c("10001_1", "10001_2", 
               "10004_1", "10004_2", 
               "10005_1", "10005_2",
               "10006_1", "10006_2",
               "10009_1", "10009_2",
               "10011_1", "10011_2",
               "10012_1", "10012_2",
               "10013_1", "10013_2",
               "10015_1", "10015_2",
               "10020_1", "10020_2",
               "10022_1", "10022_2",
               "10023_1", "10023_2",
               "10024_1", "10024_2")
Mval_all122<-Mval_all122[,col_order]

# Arrange G954_all by ID.Visit (to match column order of Mval_all122)
G954_all<-G954_all %>% arrange(ID.Visit)

# Create binary variable for Visit
G954_all$Visitdes<-ifelse(G954_all$Visit==1,0,1)
G954_all$Visitdes<-factor(G954_all$Visit,levels=c(1,2))
# Create design matrix (needed for differential methylation analysis)
design <- model.matrix(~0+G954_all$Visitdes
                       #  +celltype$Bcell
                       #  +celltype$CD4T
                       #   +celltype$CD8T
                       #  +celltype$Mono
                       #   +celltype$Neu
                       #   +celltype$NK
)
colnames(design) <- c("Visit.1", "Visit.2")
design

# Make contrast matrix since comparisons are not orthogonal (they're paired)
# contrast_matrix <- makeContrasts(Visit.2-Visit.1, levels=design)

# Test for differential methylation using the lmFit and eBayes 
# functions from limma (input data = matrix of M-values)
fit.reduced <- lmFit(Mval_all122,design)
contrast_fit<-contrasts.fit(fit.reduced,contrasts=contrast_matrix)
contrast_ebayes <- eBayes(contrast_fit, robust=TRUE)
# Extract the top differentially methylated CpGs at Visit 2
top<-topTable(contrast_ebayes, num=Inf) # negative logFC = hypomethylation in 
# Visit 2 compared to Visit 1

# Main Transmasculine Sex/Age DMP from Shepherd et al. (2022) paper
top[rownames(top)%in%"cg23256579",]
# Shepherd et al. (2022) saw significant hypomethylation at this locus after
# 12 months of transmasculinizing GAH


######## RUVm for Visit 1 and 2 ##########

# Set up the factor of interest
grp <- factor(G954_all$Visitdes)
# Extract Illumina negative control data from RGChannelSetExtended
INCs <- getINCs(rg)
# Select columns for samples of interest
INCs<-INCs[,(colnames(INCs) %in% string_use)] 
# Rename INCs columns with ID.Visit matches
colnames(INCs)<-G954_all$ID.Visit[match(colnames(INCs), 
                                               G954_all$Sample_Name)]

# Place INC columns in numerical/Visit order
col_order <- c("10001_1", "10001_2", 
               "10004_1", "10004_2", 
               "10005_1", "10005_2",
               "10006_1", "10006_2",
               "10009_1", "10009_2",
               "10011_1", "10011_2",
               "10012_1", "10012_2",
               "10013_1", "10013_2",
               "10015_1", "10015_2",
               "10020_1", "10020_2",
               "10022_1", "10022_2",
               "10023_1", "10023_2",
               "10024_1", "10024_2")
INCs<-INCs[,col_order]
#saveRDS(INCs,file="INCs.RDS")
# Add negative control data to M-values
Mc <- rbind(Mval_all122,INCs)
# Create vector marking negative controls in data matrix
ctl1 <- rownames(Mc) %in% rownames(INCs)
table(ctl1)
# Stage 1 analysis: perform RUV adjustment and fit to rank the CpGs with respect 
# to association with factor of interest
rfit1 <- RUVfit(Y = Mc, X = grp, ctl = ctl1, randomization=T)
rfit2 <- RUVadj(Y = Mc, fit = rfit1)
# Designate the CpGs that are LEAST associated with the factor of 
# interest (based on FDR-adjusted p-value) as empirical control probes (ECPs)
top1 <- topRUV(rfit2, num=Inf, p.BH = 1)
ctl2 <- rownames(Mval_all122) %in% 
    rownames(top1[top1$p.BH_X1.2 > 0.5,])
table(ctl2)
# Stage 2 analysis: use the ECPs to perform a second differential methylation 
# with RUV-inverse, which is adjusted for the unwanted variation estimated from 
# the data; perform RUV adjustment and fit
rfit3 <- RUVfit(Y = Mval_all122, X = grp, ctl = ctl2)
rfit4 <- RUVadj(Y = Mval_all122, fit = rfit3)
# Visualize top-ranked differentially methylated CpGs obtained from RUV 
# differential methylation analysis
    # p.BH = cutoff value for Benjamini-Hochberg adjusted p-values
    # p.BH_X1.2 = BH-adjusted p-values associated with factor of interest
fitdf<-as.data.frame(topRUV(rfit4,number=Inf))
saveRDS(fitdf,file="fitdf.RDS")
# Count p-values < 0.05
table(fitdf$p.BH_X1.2<0.05)
# no CpGs are significantly differentially methylated after RUVm

###### Manhattan Plot Visit 1 and 2 #########

# Annotate top
# Annotate chromosome number
top$CHR<-man$seqnames[match(rownames(top), man$Probe_ID)]
top$CHR <- str_remove(top$CHR, "chr") # Remove "chr" string
top$CHR<-as.numeric(top$CHR) # Make numeric
# Annotate chromosome position
top$BP<-man$start[match(rownames(top), man$Probe_ID)]
top$BP<-as.numeric(top$BP) # Make numeric
# Add "P" column for  P.Value
top$P<-top$P.Value
# Annotate SNPs (should be 0 because removed SNPs)
top$SNP<-man$MASK_snp5_GMAF1p[match(rownames(top), man$Probe_ID)]
top$SNP<-as.numeric(top$SNP) # Make numeric

# Make Manhattan plot using (unadjusted) p values
pdf("Manhattan 3-5 Months vs Baseline All CpG Sites.pdf")
manhattan(top, chr = "CHR", bp = "BP", p = "P", snp="SNP",
          # set threshold as Bonferroni p-value ??
          suggestiveline = -log10(2.5e-6), #blue; determined arbitrarily
          genomewideline= -log10(0.05/length(top$P.Value)), #red; Bonferroni
          # colors
          col = c("gray66", "gray45"),
          # annotateTop=TRUE,
          main="3-5 Months vs Baseline All CpG Sites",
          cex=0.7,
          ylim=c(0,8)) 
dev.off()

topsig<-top[top$P.Value<0.05,]

############# Visit 1 and 2 on Testosterone #############

# Create column in metadata that combines Subject.ID with Visit
G954_allTB<-unite(G954_allTB, col="ID.Visit", c("Subject.ID","Visit"), sep="_",
                remove=F)
# Rename Mval_all122 columns with ID.Visit matches
colnames(Mval_all122TB)<-G954_allTB$ID.Visit[match(colnames(Mval_all122TB), 
                                               G954_allTB$Sample_Name)]

# Place Mval columns in numerical/Visit order
col_orderTB <- c("10004_1", "10004_2", 
               "10005_1", "10005_2",
               "10006_1", "10006_2",
               "10009_1", "10009_2",
               "10011_1", "10011_2",
               "10012_1", "10012_2",
               "10015_1", "10015_2",
               "10020_1", "10020_2",
               "10022_1", "10022_2",
               "10024_1", "10024_2")
Mval_all122TB<-Mval_all122TB[,col_orderTB]

# Arrange G954_all by ID.Visit (to match column order of Mval_all122)
G954_allTB<-G954_allTB %>% arrange(ID.Visit)

# Create binary variable for Visit
G954_allTB$Visitdes<-ifelse(G954_allTB$Visit==1,0,1)
G954_allTB$Visitdes<-factor(G954_allTB$Visit,levels=c(1,2))
# Create design matrix (needed for differential methylation analysis)
designTB <- model.matrix(~0+G954_allTB$Visitdes
                       #  +celltype$Bcell
                       #  +celltype$CD4T
                       #   +celltype$CD8T
                       #  +celltype$Mono
                       #   +celltype$Neu
                       #   +celltype$NK
)
colnames(designTB) <- c("Visit.1", "Visit.2")
designTB

# Make contrast matrix since comparisons are not orthogonal (they're paired)
contrast_matrixTB <- makeContrasts(Visit.2-Visit.1, levels=designTB)

# Test for differential methylation using the lmFit and eBayes 
# functions from limma (input data = matrix of M-values)
fit.reducedTB <- lmFit(Mval_all122TB,designTB)
contrast_fitTB<-contrasts.fit(fit.reducedTB,contrasts=contrast_matrixTB)
contrast_ebayesTB <- eBayes(contrast_fitTB, robust=TRUE)
# Extract the top differentially methylated CpGs at Visit 2
topTB<-topTable(contrast_ebayesTB, num=Inf) # negative logFC = hypomethylation in 
# Visit 2 compared to Visit 1
# adj p value = 1, even when leaving X and Y in

# Annotate topTB
# Annotate chromosome number
topTB$CHR<-man$seqnames[match(rownames(topTB), man$Probe_ID)]
topTB$CHR <- str_remove(topTB$CHR, "chr") # Remove "chr" string
topTB$CHR<-as.numeric(topTB$CHR) # Make numeric
# Annotate chromosome position
topTB$BP<-man$start[match(rownames(topTB), man$Probe_ID)]
topTB$BP<-as.numeric(topTB$BP) # Make numeric
# Add "P" column for  P.Value
topTB$P<-topTB$P.Value
# Annotate SNPs (should be 0 because removed SNPs)
topTB$SNP<-man$MASK_snp5_GMAF1p[match(rownames(topTB), man$Probe_ID)]
topTB$SNP<-as.numeric(topTB$SNP) # Make numeric

# Make Manhattan plot using (unadjusted) p values
pdf("Manhattan 3-5 Months vs Baseline All CpG Sites (Transmasculine).pdf")
manhattan(topTB, chr = "CHR", bp = "BP", p = "P", snp="SNP",
          # set threshold as Bonferroni p-value ??
          suggestiveline = -log10(2.5e-6), #blue; how to determine this?
          genomewideline= -log10(0.05/length(topTB$P.Value)), #red; Bonferroni
          # colors
          col = c("gray66", "gray45"),
          # annotateTop=TRUE,
          main="3-5 Months vs Baseline All CpG Sites (Transmasculine Youth)",
          cex=0.7,
          ylim=c(0,8)) 
dev.off()

topTBsig<-topTB[topTB$P.Value<0.05,]

########### Visit 1 and 3 #############
# Create column in metadata that combines Subject.ID with Visit
G954_13<-unite(G954_13, col="ID.Visit", c("Subject.ID","Visit"), sep="_",
                remove=F)
# Rename Mval_all122 columns with ID.Visit matches
colnames(Mval_all133)<-G954_13$ID.Visit[match(colnames(Mval_all133), 
                                               G954_13$Sample_Name)]

# Place Mval columns in numerical/Visit order
col_order13 <- c("10001_1", "10001_3", 
               "10004_1", "10004_3", 
               "10005_1", "10005_3")
Mval_all133<-Mval_all133[,col_order13]

# Arrange G954_all by ID.Visit (to match column order of Mval_all133)
G954_13<-G954_13 %>% arrange(ID.Visit)

# Create binary variable for Visit
G954_13$Visitdes<-ifelse(G954_13$Visit==1,0,1)
G954_13$Visitdes<-factor(G954_13$Visit,levels=c(1,3))
# Create design matrix (needed for differential methylation analysis)
design13 <- model.matrix(~0+G954_13$Visitdes
                       #  +celltype$Bcell
                       #  +celltype$CD4T
                       #   +celltype$CD8T
                       #  +celltype$Mono
                       #   +celltype$Neu
                       #   +celltype$NK
)
colnames(design13) <- c("Visit.1", "Visit.3")
design13

# Make contrast matrix since comparisons are not orthogonal (they're paired)
contrast_matrix13 <- makeContrasts(Visit.3-Visit.1, levels=design13)

# Test for differential methylation using the lmFit and eBayes 
# functions from limma (input data = matrix of M-values)
fit.reduced13 <- lmFit(Mval_all133,design13)
contrast_fit13<-contrasts.fit(fit.reduced13,contrasts=contrast_matrix13)
contrast_ebayes13 <- eBayes(contrast_fit13, robust=TRUE)
# Extract the top differentially methylated CpGs at Visit 3
top13<-topTable(contrast_ebayes13, num=Inf) # negative logFC = hypomethylation in 
# Visit 3 compared to Visit 1

# Annotate top
# Annotate chromosome number
top13$CHR<-man$seqnames[match(rownames(top13), man$Probe_ID)]
top13$CHR <- str_remove(top13$CHR, "chr") # Remove "chr" string
top13$CHR<-as.numeric(top13$CHR) # Make numeric
# Annotate chromosome position
top13$BP<-man$start[match(rownames(top13), man$Probe_ID)]
top13$BP<-as.numeric(top13$BP) # Make numeric
# Add "P" column for P.Value
top13$P<-top13$P.Value
# Annotate SNPs (should be 0 because removed SNPs)
top13$SNP<-man$MASK_snp5_GMAF1p[match(rownames(top13), man$Probe_ID)]
top13$SNP<-as.numeric(top13$SNP) # Make numeric

# Make Manhattan plot using (unadjusted) p values
pdf("Manhattan 12 Months vs Baseline All CpG Sites.pdf")
manhattan(top13, chr = "CHR", bp = "BP", p = "P", snp="SNP",
          # set threshold as Bonferroni p-value (0.05/# CpGs)
          suggestiveline = -log10(2.5e-6), 
          genomewideline= -log10(0.05/length(top13$P.Value)),
          # colors
          col = c("gray66", "gray45"),
          # annotateTop=TRUE,
          main="12 Months vs Baseline All CpG Sites",
          cex=0.7,
          ylim=c(0,8)) 
dev.off()

top13sig<-top13[top13$P.Value<0.05,]


########### Visit 1 and 3 on Testosterone #############
# Create column in metadata that combines Subject.ID with Visit
G954_13TB<-unite(G954_13TB, col="ID.Visit", c("Subject.ID","Visit"), sep="_",
               remove=F)
# Rename Mval_all122 columns with ID.Visit matches
colnames(Mval_all133TB)<-G954_13TB$ID.Visit[match(colnames(Mval_all133TB), 
                                              G954_13TB$Sample_Name)]

# Place Mval columns in numerical/Visit order
col_order13TB <- c("10004_1", "10004_3", 
                 "10005_1", "10005_3")
Mval_all133TB<-Mval_all133TB[,col_order13TB]

# Arrange G954_all by ID.Visit (to match column order of Mval_all133)
G954_13TB<-G954_13TB %>% arrange(ID.Visit)

# Create binary variable for Visit
G954_13TB$Visitdes<-ifelse(G954_13TB$Visit==1,0,1)
G954_13TB$Visitdes<-factor(G954_13TB$Visit,levels=c(1,3))
# Create design matrix (needed for differential methylation analysis)
design13TB <- model.matrix(~0+G954_13TB$Visitdes
                         #  +celltype$Bcell
                         #  +celltype$CD4T
                         #   +celltype$CD8T
                         #  +celltype$Mono
                         #   +celltype$Neu
                         #   +celltype$NK
)
colnames(design13TB) <- c("Visit.1", "Visit.3")
design13TB

# Make contrast matrix since comparisons are not orthogonal (they're paired)
contrast_matrix13TB <- makeContrasts(Visit.3-Visit.1, levels=design13TB)

# Test for differential methylation using the lmFit and eBayes 
# functions from limma (input data = matrix of M-values)
fit.reduced13TB <- lmFit(Mval_all133TB,design13TB)
contrast_fit13TB<-contrasts.fit(fit.reduced13TB,contrasts=contrast_matrix13TB)
contrast_ebayes13TB <- eBayes(contrast_fit13TB, robust=TRUE)
# Extract the top differentially methylated CpGs at Visit 3
top13TB<-topTable(contrast_ebayes13TB, num=Inf) # negative logFC = hypomethylation in 
# Visit 3 compared to Visit 1

# Annotate top13TB
# Annotate chromosome number
top13TB$CHR<-man$seqnames[match(rownames(top13TB), man$Probe_ID)]
top13TB$CHR <- str_remove(top13TB$CHR, "chr") # Remove "chr" string
top13TB$CHR<-as.numeric(top13TB$CHR) # Make numeric
# Annotate chromosome position
top13TB$BP<-man$start[match(rownames(top13TB), man$Probe_ID)]
top13TB$BP<-as.numeric(top13TB$BP) # Make numeric
# Add "P" column for P.Value
top13TB$P<-top13TB$P.Value
# Annotate SNPs (should be 0 because removed SNPs)
top13TB$SNP<-man$MASK_snp5_GMAF1p[match(rownames(top13TB), man$Probe_ID)]
top13TB$SNP<-as.numeric(top13TB$SNP) # Make numeric

# Make Manhattan plot using (unadjusted) p values
pdf("Manhattan 12 Months vs Baseline All CpG Sites (Transmasculine).pdf")
manhattan(top13TB, chr = "CHR", bp = "BP", p = "P", snp="SNP",
          # set threshold as Bonferroni p-value (0.05/# CpGs)
          suggestiveline = -log10(2.5e-6), 
          genomewideline= -log10(0.05/length(top13TB$P.Value)),
          # colors
          col = c("gray66", "gray45"),
          # annotateTop=TRUE,
          main="12 Months vs Baseline All CpG Sites (Transmasculine Youth)",
          cex=0.7,
          ylim=c(0,8)) 
dev.off()

top13TBsig<-top13TB[top13TB$P.Value<0.05,]

# Check with top differentially methylated CpG from Shepherd et al. (2022)
top13TB[rownames(top13TB)%in%"cg23256579",]
# Shepherd et al. (2022) saw significant hypomethylation at this locus after
# 12 months of transmasculinizing GAH


###### Gene Annotation #####
# Pick out the CpGs above the blue suggestiveline
top13TBsuggestive<-top13TBsig[(top13TBsig$P.Value < 2.5e-6),]
# Pick out these CpGs in Visit 1/2 TB data
topTBsuggestive<-topTB[(rownames(topTB)%in%rownames(top13TBsuggestive)),]

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets")
# Load Illumina manifest file
illumina.man<-read.csv("MethylationEPIC_v-1-0_B4.csv")

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/EWAS")

# Get Illumina manifest info for the top 2 CpGs
sigCpGs_illumina<-illumina.man[(illumina.man$IlmnID %in% c("cg14683750",
                                             "cg16655805",
                                             "cg20386316",
                                             "cg06916709",
                                             "cg05126421",
                                             "cg21895450")),]

# View relevant annotation columns
sigCpGs_illumina[,colnames(sigCpGs_illumina) %in% 
                     c("IlmnID", "CHR",
                       "UCSC_RefGene_Name", 
                       "UCSC_RefGene_Group", 
                       "UCSC_CpG_Islands_Name",
                       "Relation_to_UCSC_CpG_Island",
                       "GencodeBasicV12_NAME", 
                       "GencodeBasicV12_Group", 
                       "GencodeCompV12_NAME", 
                       "GencodeCompV12_Group",
                       "DMR",
                       "Regulatory_Feature_Group")] %>% View()

# Get manifest info for the top 6 CpGs
sigCpGs_suggestive<-man[(man$Probe_ID %in% c("cg14683750",
                         "cg16655805",
                         "cg20386316",
                         "cg06916709",
                         "cg05126421",
                         "cg21895450")),] %>% select(seqnames, start, end, 
                                                strand, gene, gene_HGNC)

# Get manifest for CpGs annotated to PRR4 gene (from Shepherd et al. 2022)
man[(man$gene_HGNC %in% "PRR4"),] # no CpGs are annotated to PRR4 gene

##### Comparisons with Shepherd et al. (2022) findings #####
# Transmasculine Sex/Age DMPs from Shepherd et al. (2022) paper
TBshepherd<-topTB[rownames(topTB)%in%c("cg07658646",
                           "cg23256579", # Most significant Shepherd finding
                           "cg27615582",
                           "cg13942826"),]
# Shepherd et al. (2022) saw significant hypomethylation at these loci after
# 12 months of transmasculinizing GAH

# Transmasculine Sex/Age DMPs from Shepherd et al. (2022) paper
TB13shepherd<-top13TB[rownames(top13TB)%in%c("cg07658646",
                               "cg23256579", # Most significant Shepherd finding
                               "cg27615582",
                               "cg13942826"),]
# Shepherd et al. (2022) saw significant hypomethylation at these loci after
# 12 months of transmasculinizing GAH

# Match the order of rows using row names
topshepherd_matchedTB <- TBshepherd[match(row.names(TB13shepherd), 
                                                 row.names(TBshepherd)), ]

# Differentiate logFC columns
topshepherd_matchedTB$logFC_12<-topshepherd_matchedTB$logFC
TB13shepherd$logFC_13<-TB13shepherd$logFC
# Differentiate P.Value columns
topshepherd_matchedTB$P.Value_12<-topshepherd_matchedTB$P.Value
TB13shepherd$P.Value_13<-TB13shepherd$P.Value

# Combine Visit 1/2 and Visit 1/3 significant (unadjusted) CpGs
top_shepherdTB<-cbind(topshepherd_matchedTB,TB13shepherd)

##### Make plot in same style of Shepherd (2022) plot #####
# Convert M-values to beta values
beta_all122TB<-M2B(Mval_all122TB)
beta_all133TB<-M2B(Mval_all133TB)

# Separate by time point
beta_baselineTB<-beta_all122TB[,str_detect(colnames(beta_all122TB),"_1")]
beta_2TB<-beta_all122TB[,str_detect(colnames(beta_all122TB),"_2")]
beta_3TB<-beta_all133TB[,str_detect(colnames(beta_all133TB),"_3")]

# Match row order to baseline row order using row names
beta_2TB <- beta_2TB[match(row.names(beta_baselineTB), 
                                         row.names(beta_2TB)), ]
beta_3TB <- beta_3TB[match(row.names(beta_baselineTB), 
                           row.names(beta_3TB)), ]

# Select transmasculine DMPs from Shepherd et al. (2022) paper
beta_baselineTBshep<-beta_baselineTB[rownames(beta_baselineTB)%in%c("cg07658646",
                                             "cg23256579", # Most significant Shepherd finding
                                             "cg27615582",
                                             "cg13942826"),]
beta_2TBshep<-beta_2TB[rownames(beta_2TB)%in%c("cg07658646",
                                               "cg23256579", # Most significant Shepherd finding
                                               "cg27615582",
                                               "cg13942826"),]
beta_3TBshep<-beta_3TB[rownames(beta_3TB)%in%c("cg07658646",
                                               "cg23256579", # Most significant Shepherd finding
                                               "cg27615582",
                                               "cg13942826"),]

# Rearrange rownames to match order in Shepherd (2022) plot
# Define new order of row names
new_order <- c("cg07658646",
               "cg23256579", # Most significant Shepherd finding
               "cg27615582",
               "cg13942826")

# Reorder row names
beta_baselineTBshep <- beta_baselineTBshep[new_order, ]
beta_2TBshep <- beta_2TBshep[new_order, ]
beta_3TBshep <- beta_3TBshep[new_order, ]

# Make data frame
beta_baselineTBshep <- as.data.frame((beta_baselineTBshep))
beta_2TBshep <- as.data.frame((beta_2TBshep))
beta_3TBshep <- as.data.frame((beta_3TBshep))

# Transpose 
beta_baset <- t(beta_baselineTBshep)
beta_2t <- t(beta_2TBshep)
beta_3t <- t(beta_3TBshep)

# Make dataframe
beta_baset <- as.data.frame(beta_baset)
beta_2t <- as.data.frame(beta_2t)
beta_3t <- as.data.frame(beta_3t)

# Beta Value Statistics
# Perform Shapiro test to see distribution
shapiro.test(x=c(beta_baset$cg07658646,
                 beta_2t$cg07658646)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(beta_baset$cg07658646,
       beta_2t$cg07658646,paired=T)
# Perform Shapiro test to see distribution
shapiro.test(x=c(beta_baset$cg23256579,
                 beta_2t$cg23256579)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(beta_baset$cg23256579,
       beta_2t$cg23256579,paired=T)
# Perform Shapiro test to see distribution
shapiro.test(x=c(beta_baset$cg27615582,
                 beta_2t$cg27615582)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(beta_baset$cg27615582,
       beta_2t$cg27615582,paired=T)
# Perform Shapiro test to see distribution
shapiro.test(x=c(beta_baset$cg13942826,
                 beta_2t$cg13942826)) 
# Perform paired t-test (or Wilcoxon in Prism if distribution not normal)
t.test(beta_baset$cg13942826,
       beta_2t$cg13942826,paired=T)

# Perform Shapiro test to see distribution
shapiro.test(x=c(beta_2t$cg07658646,
                 beta_3t$cg07658646)) 
# Perform paired t-test (or Mann-Whitney in Prism if distribution not normal)
t.test(beta_2t$cg07658646,
       beta_3t$cg07658646,paired=F)
# Perform Shapiro test to see distribution
shapiro.test(x=c(beta_2t$cg23256579,
                 beta_3t$cg23256579)) 
# Perform paired t-test (or Mann-Whitney in Prism if distribution not normal)
t.test(beta_2t$cg23256579,
       beta_3t$cg23256579,paired=F)
# Perform Shapiro test to see distribution
shapiro.test(x=c(beta_2t$cg27615582,
                 beta_3t$cg27615582)) 
# Perform paired t-test (or Mann-Whitney in Prism if distribution not normal)
t.test(beta_2t$cg27615582,
       beta_3t$cg27615582,paired=F)
# Perform Shapiro test to see distribution
shapiro.test(x=c(beta_2t$cg13942826,
                 beta_3t$cg13942826)) 
# Perform paired t-test (or Mann-Whitney in Prism if distribution not normal)
t.test(beta_2t$cg13942826,
       beta_3t$cg13942826,paired=F)


# Make new columns with row means and standard deviations
beta_baselineTBshep$mean_beta1 <- rowMeans(beta_baselineTBshep)
beta_baselineTBshep$sd_beta1 <- apply(beta_baselineTBshep,1,sd)
beta_2TBshep$mean_beta2 <- rowMeans(beta_2TBshep)
beta_2TBshep$sd_beta2 <- apply(beta_2TBshep,1,sd)
beta_3TBshep$mean_beta3 <- rowMeans(beta_3TBshep)
beta_3TBshep$sd_beta3 <- apply(beta_3TBshep,1,sd)

# Combine Visit 1,2, and 3 beta data
beta_allTBshep<-cbind(beta_baselineTBshep,beta_2TBshep,beta_3TBshep)
# Select only the mean and sd columns
beta_allTBshep_meansd<-beta_allTBshep[,colnames(beta_allTBshep) %in%
                                        c("mean_beta1",
                                          "sd_beta1",
                                          "mean_beta2",
                                          "sd_beta2",
                                          "mean_beta3",
                                          "sd_beta3")]

# Save as CSV
write.csv(beta_allTBshep_meansd,file="TBshep_meansd.csv")