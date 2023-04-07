# Differential Methylation Analysis with RUVm

# Load packages
library(dplyr)
library(missMethyl)
library(ENmix)
library(limma)
library(ggplot2)
library(DMRcate)

# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

###### DMP and DMR Identification in Transmasculine Youth ######
# Load objects from "pre-processing and normalization.R" code
# G954 metadata
G954sheet <- readRDS("G954sheet.RDS")
# G954 RGChannelSetExtended
rg <- readRDS("G954_EPIC_Data_0217rgsetext.RDS")
# G954 beta values
rg_betas <- readRDS("G954_EPIC_Data_0217rg_betas.RDS")

# Convert beta values to M-values using ENmix
MvalAFAB<-B2M(rg_betas)

# Remove columns for people not of interest (i.e., unusable paired data,
# people only on blockers, data from Visit 3, or people assigned male at birth)
G954_GAH<-G954sheet %>% 
    subset(., !(Sample_Name %in% c("206839580019_R01C01",# 10002
                                   "206839580019_R02C01", # 10002_2
                                   "206839580032_R07C01",# 10017
                                   "206839580032_R08C01",# 10017_2
                                   "206839580066_R02C01",# 10018
                                   "206839580066_R01C01",# 10018_2
                                   "206839580077_R08C01",#10007_2
                                   "206839580077_R07C01",#10007
                                   "206839580077_R03C01",#10005_3
                                   "206839580077_R04C01",#10004_3
                                   "206842050024_R02C01", #10001_3 
                                   "206842050024_R03C01", #10001
                                   "206842050024_R01C01", #10001_2
                                   "206839580032_R03C01", #10013
                                   "206839580032_R04C01", #10013_2
                                   "206839580066_R05C01", #10023
                                   "206839580066_R06C01")))  #10023_2

# Create design matrix
groupAFAB <- factor(G954_GAH$Visit,levels=c(1,2))
idAFAB <- factor(G954_GAH$Subject.ID)
designAFAB <- model.matrix(~idAFAB + groupAFAB)
designAFAB

# Create character vector of Sample_Names of interest
AFAB_string_use<-G954_GAH %>% pull(Sample_Name)

# Select columns for samples of interest (people assigned female at birth)
MvalAFAB_no3<-MvalAFAB[,(colnames(MvalAFAB) %in% AFAB_string_use)] 
saveRDS(MvalAFAB_no3,file="MvalAFAB_no3.RDS")

# Test for differential methylation using the lmFit and eBayes 
# functions from limma (input data = matrix of M-values)
fit.reducedAFAB <- lmFit(MvalAFAB_no3,designAFAB)
fit.reducedAFAB <- eBayes(fit.reducedAFAB, robust=TRUE)
# Check the numbers of hyper-methylated (1) and hypo-methylated (-1) 
# using the decideTests function in limma 
summaryAFAB<-summary(decideTests(fit.reducedAFAB))
saveRDS(summaryAFAB,file="summaryAFAB.RDS")
decideTestsAFAB<-decideTests(fit.reducedAFAB)
saveRDS(decideTestsAFAB,file="decideTestsAFAB.RDS")
# Extract the top 10 differentially methylated CpGs for Visit 1 versus Visit 2 
# using topTable
topAFAB<-topTable(fit.reducedAFAB,coef=11)
topAFAB
saveRDS(topAFAB,file="topAFAB.RDS")
# Set up the factor of interest
grpAFAB <- factor(G954_GAH$Visit, labels=c(1,2))
# Extract Illumina negative control data from RGChannelSetExtended
INCsAFAB <- getINCs(rg)
# Select columns for samples of interest
INCsAFAB_no3<-INCsAFAB[,(colnames(INCsAFAB) %in% AFAB_string_use)] 
saveRDS(INCsAFAB_no3,file="INCsAFAB_no3.RDS")
# Add negative control data to M-values
McAFAB <- rbind(MvalAFAB_no3,INCsAFAB_no3)
# Create vector marking negative controls in data matrix
ctl1AFAB <- rownames(McAFAB) %in% rownames(INCsAFAB_no3)
table(ctl1AFAB)
# Stage 1 analysis: perform RUV adjustment and fit to rank the CpGs with respect 
# to association with factor of interest
rfit1AFAB <- RUVfit(Y = McAFAB, X = grpAFAB, ctl = ctl1AFAB)
rfit2AFAB <- RUVadj(Y = McAFAB, fit = rfit1AFAB)
# Designate the CpGs that are LEAST associated with the factor of 
# interest (based on FDR-adjusted p-value) as empirical control probes (ECPs)
top1AFAB <- topRUV(rfit2AFAB, num=Inf, p.BH = 1)
ctl2AFAB <- rownames(MvalAFAB_no3) %in% 
    rownames(top1AFAB[top1AFAB$p.BH_X1.2 > 0.5,])
table(ctl2AFAB)
# Stage 2 analysis: use the ECPs to perform a second differential methylation 
# with RUV-inverse, which is adjusted for the unwanted variation estimated from 
# the data; perform RUV adjustment and fit
rfit3AFAB <- RUVfit(Y = MvalAFAB_no3, X = grpAFAB, ctl = ctl2AFAB)
rfit4AFAB <- RUVadj(Y = MvalAFAB_no3, fit = rfit3AFAB)
# Visualize top-ranked differentially methylated CpGs obtained from RUV 
# differential methylation analysis
    # p.BH = cutoff value for Benjamini-Hochberg adjusted p-values
    # p.BH_X1.2 = BH-adjusted p-values associated with factor of interest
fitdfAFAB<-as.data.frame(topRUV(rfit4AFAB,number=Inf))
saveRDS(fitdfAFAB,file="fitdfAFAB.RDS")
# Count p values < 0.05
table(fitdfAFAB$p.BH_X1.2<0.05)

# Select all significant CpGs in dataset
top9357<-as.data.frame(topRUV(rfit4AFAB,number=9357))
saveRDS(top9357,file="top9357.RDS")
top9357<-readRDS("top9357.RDS")
sigCpGsAFAB<-rownames(top9357)
# Check number of genes to which significant CpGs are annotated
checkAFAB <- getMappedEntrezIDs(sig.cpg = sigCpGsAFAB)
length(checkAFAB$sig.eg)
# ^ 9357 CpGs are differentially methylated in GAH participants assigned female
# at birth between Visits 1 and 2, mapping to 2572 different genes

# Perform gene set enrichment analysis of DMPs with missMethyl
gstAFAB <- gometh(sig.cpg=sigCpGsAFAB, all.cpg=rownames(fitdfAFAB), 
                  collection="GO",plot.bias=TRUE, array.type="EPIC",sig.genes=F)
topGSA.gstAFAB<-topGSA(gstAFAB, n=Inf)
saveRDS(topGSA.gstAFAB,file="topGSA.gst.AFAB.4-5-23.RDS")

# Create a bubble plot of the top 10 GO terms for DMPs
topGSA.gstAFAB<-readRDS("topGSA.gst.AFAB.4-5-23.RDS")
top10AFAB<-topGSA.gstAFAB %>% head(n=10)
top10AFAB %>% select(TERM)
ggplot(top10AFAB, aes(x = DE, y = TERM, size = -log10(P.DE), label =(""))) +
    geom_point(alpha = 0.7, color = "red") +
    geom_text(size = 5, vjust = 2) +
    scale_size(range = c(2, 10)) +
    labs(x = "Differentially Methylated Genes", y = "", 
         title = "Top 10 GO Terms for DMPs in Transmasculine Youth on GAH",
         size = "-log10(P-value)") +
    theme_bw() + 
    theme(axis.text = element_text(size = 17), 
          axis.title.x=element_text(size = 20,
                                    margin = margin(t = 30, r = 0, b = 10, 
                                                    l = 0), 
                                    vjust = 2),
          legend.title=element_text(size = 16),
          legend.text=element_text(size = 14))+
    scale_x_continuous(limits = c(0, 19000))

# Annotate M-values with genomic information in DMRcate
myAnnotationAFAB <- cpg.annotate(object = MvalAFAB_no3, datatype = "array", 
                                 what = "M", 
                                 arraytype = c("EPIC"), 
                                 analysis.type = "differential", 
                                 design = designAFAB, 
                                 coef = 11)
# Identify DMRs; create DMResults object
DMRsAFAB <- dmrcate(myAnnotationAFAB, lambda=1000, C=2) # 4 DMRs identified
# Extract DMR range information; create GRanges object
results.rangesAFAB <- extractRanges(DMRsAFAB)
results.rangesAFAB
saveRDS(results.rangesAFAB, file="results.rangesAFAB.RDS")
write.csv(results.rangesAFAB,file="results.rangesAFAB.csv")

# Perform gene set enrichment analysis of DMRs with missMethyl
gst.regionAFAB <- goregion(results.rangesAFAB, all.cpg=rownames(MvalAFAB_no3), 
                           collection="GO", array.type="EPIC", plot.bias=TRUE,
                           sig.genes=F
)
topGSA.gst.regionAFAB<-topGSA(gst.regionAFAB, n=Inf)
saveRDS(topGSA.gst.regionAFAB,file="topGSA.gst.region.AFAB.4-5-23.RDS")

# Create a bubble plot of the top 10 GO terms for DMRs
topGSA.gst.regionAFAB<-readRDS("topGSA.gst.region.AFAB.4-5-23.RDS")
top10.gst.regionAFAB<-topGSA.gst.regionAFAB %>% head(n=10)
ggplot(top10.gst.regionAFAB, aes(x = DE, y = TERM, size = -log10(P.DE), 
                                 label =(""))) +
    geom_point(alpha = 0.7, color = "red") +
    geom_text(size = 5, vjust = 2) +
    scale_size(range = c(2, 10)) +
    labs(x = "Differentially Methylated Genes", y = "", 
         title = "Top 10 GO Terms for DMRs in Transmasculine Youth on GAH",
         size = "-log10(P-value)") +
    theme_bw() +
    theme(axis.text = element_text(size = 20), 
          axis.title.x=element_text(size = 20,
                                    margin = margin(t = 30, r = 0, b = 10, 
                                                    l = 0), 
                                    vjust = 2),
          axis.title.y=element_text(size = 20,
                                    margin = margin(t = 0, r = 75, b = 0, 
                                                    l = 0), 
                                    vjust = 2),
          legend.title=element_text(size = 16),
          legend.text=element_text(size = 14))+
    scale_x_continuous(limits = c(0, 19000))


###### DMP and DMR Identification in Transmasculine Adults (GSE173717) ######
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/03-28-23")

# Load objects from "pre-processing and normalization.R" code
# GSE173717 metadata
sheet <- readRDS("GSE173717_Sample_Sheet.RDS")
# GSE173717 RGChannelSetExtended
rgGSE <- readRDS("GSE173717rgsetext.RDS")
# GSE173717 beta values
rg_betasGSE <- readRDS("GSE173717rg_betas.RDS")

# Convert beta values to M-values using ENmix
MvalGSE<-B2M(rg_betasGSE)

# Select columns for subjects/visits of interest (transmasculine, Visits 1 & 2)
GSE173717_GAH<-sheet %>% 
    subset(., (Visit==1&Female==1&Gender.Identity=="Transgender"|
                   Visit==2&Female==1&Gender.Identity=="Transgender"))

# Create character vector of Sample_Names not of interest
GSE_string_dontuse<-sheet %>% 
    subset(., !(Visit==1&Female==1&Gender.Identity=="Transgender"|
                    Visit==2&Female==1&Gender.Identity=="Transgender")) %>% 
    pull(Sample_Name)

# Create design matrix
groupGSE <- factor(GSE173717_GAH$Visit,levels=c(1,2))
saveRDS(groupGSE,file="groupGSE.RDS")
idGSE <- factor(GSE173717_GAH$Subject.ID)
designGSE <- model.matrix(~idGSE + groupGSE)
designGSE

# Remove columns for samples not of interest
MvalGSE_no3<-MvalGSE[,!(colnames(MvalGSE) %in% GSE_string_dontuse)] 
saveRDS(MvalGSE_no3,file="MvalGSE_no3.RDS")
# Test for differential methylation using the lmFit and eBayes 
# functions from limma (input data = matrix of M-values)
fit.reducedGSE <- lmFit(MvalGSE_no3,designGSE)
fit.reducedGSE <- eBayes(fit.reducedGSE, robust=TRUE)
# Check the numbers of hyper-methylated (1) and hypo-methylated (-1) 
# using the decideTests function in limma 
summary(decideTests(fit.reducedGSE))
summaryGSE<-summary(decideTests(fit.reducedGSE))
saveRDS(summaryGSE,file="summaryGSE.RDS")
# Extract the top 10 differentially methylated CpGs for Visit 1 versus Visit 2  
# using topTable
topGSE<-topTable(fit.reducedGSE,coef=10)
topGSE
saveRDS(topGSE,file="topGSE.RDS")

# Set up the factor of interest
grpGSE <- factor(GSE173717_GAH$Visit, labels=c(1,2))
# Extract Illumina negative control data
INCsGSE <- getINCs(rgGSE) 
# Remove columns for samples not of interest
INCsGSE_no3<-INCsGSE[,!(colnames(INCsGSE) %in% GSE_string_dontuse)]
# Add negative control data to M-values
McGSE <- rbind(MvalGSE_no3,INCsGSE_no3)
# Create vector marking negative controls in data matrix
ctl1GSE <- rownames(McGSE) %in% rownames(INCsGSE_no3)
table(ctl1GSE)
# Stage 1 analysis: perform RUV adjustment and fit to rank the CpGs with respect 
# to association with factor of interest
rfit1GSE <- RUVfit(Y = McGSE, X = grpGSE, ctl = ctl1GSE) # Stage 1 analysis
rfit2GSE <- RUVadj(Y = McGSE, fit = rfit1GSE)
# Designate the CpGs that are LEAST associated with the factor of 
# interest (based on FDR-adjusted p-value) as empirical control probes (ECPs)
top1GSE <- topRUV(rfit2GSE, num=Inf, p.BH = 1)
ctl2GSE <- rownames(MvalGSE_no3) %in% 
    rownames(top1GSE[top1GSE$p.BH_X1.2 > 0.5,])
table(ctl2GSE)
#We can then use the ECPs to perform a second differential methylation with 
#RUV-inverse, which is adjusted for the unwanted variation estimated from the 
#data.
# Stage 2 analysis: use the ECPs to perform a second differential methylation 
# with RUV-inverse, which is adjusted for the unwanted variation estimated from 
# the data; perform RUV adjustment and fit
rfit3GSE <- RUVfit(Y = MvalGSE_no3, X = grpGSE, ctl = ctl2GSE)# Stage 2 analysis
rfit4GSE <- RUVadj(Y = MvalGSE_no3, fit = rfit3GSE)
# Visualize top-ranked differentially methylated CpGs obtained from RUV 
# differential methylation analysis
    # p.BH = cutoff value for Benjamini-Hochberg adjusted p-values
    # p.BH_X1.2 = BH-adjusted p-values associated with factor of interest
fitdfGSE<-as.data.frame(topRUV(rfit4GSE,number=Inf))
saveRDS(fitdfGSE,file="fitdfGSE.RDS")
# Count p values < 0.05
table(fitdfGSE$p.BH_X1.2<0.05)

# Select all significant CpGs in dataset
top285255<-as.data.frame(topRUV(rfit4GSE,number=285255))
sigCpGsGSE<-rownames(top285255)
# Check number of genes to which significant CpGs are annotated
checkGSE <- getMappedEntrezIDs(sig.cpg = sigCpGsGSE)
length(checkGSE$sig.eg)
# ^ 285255 CpGs are differentially methylated in GAH participants assigned 
# female at birth between Visits 1 and 2, mapping to 19107 different genes

# Perform gene set enrichment analysis of DMPs with missMethyl
gstGSE <- gometh(sig.cpg=sigCpGsGSE, all.cpg=rownames(fitdfGSE), 
                 collection="GO",plot.bias=TRUE, array.type="EPIC",sig.genes=F)
# ^ Hypergeometric test performed for gene sets where p = exactly 0
topGSA.gst_GSE<-topGSA(gstGSE, n=Inf)
saveRDS(topGSA.gst_GSE,file="topGSA.gst.AFAB.GSE173717.4-5-23.RDS")

# Create a bubble plot of the top 10 GO terms for DMPs
topGSA.gst_GSE<-readRDS("topGSA.gst.AFAB.GSE173717.4-5-23.RDS")
top10_GSE173717<-topGSA.gst_GSE %>% head(n=10)
ggplot(top10_GSE173717, aes(x = DE, y = TERM, size = -log10(P.DE), 
                            label =(""))) +
    geom_point(alpha = 0.7, color = "blue") +
    geom_text(size = 5, vjust = 2) +
    scale_size(range = c(2, 10)) +
    labs(x = "Differentially Methylated Genes", y = "", 
         title = "Top 10 GO Terms for DMPs in Transmasculine Adults",
         size = "-log10(P-value)") +
    theme_bw() +
    theme(axis.text = element_text(size = 20), 
          axis.title.x=element_text(size = 20,
                                    margin = margin(t = 30, r = 0, b = 10, 
                                                    l = 0), 
                                    vjust = 2),
          axis.title.y=element_text(size = 20,
                                    margin = margin(t = 0, r = 250, b = 0, 
                                                    l = 0), 
                                    vjust = 2),
          legend.title=element_text(size = 16),
          legend.text=element_text(size = 14))+
    scale_x_continuous(limits = c(0, 19000))

# Annotate M-values with genomic information in DMRcate
myAnnotationGSE <- cpg.annotate(object = MvalGSE_no3, datatype = "array", 
                                what = "M", 
                                arraytype = c("EPIC"), 
                                analysis.type = "differential", 
                                design = designGSE, 
                                coef = 10)
# Identify DMRs; create DMResults object
DMRsGSE <- dmrcate(myAnnotationGSE, lambda=1000, C=2) # 59016 DMRs identified
# Extract DMR range information; create GRanges object
results.rangesGSE <- extractRanges(DMRsGSE)
sum(results.rangesGSE$Fisher<0.05)
saveRDS(results.rangesGSE, file="results.rangesGSE.RDS")
write.csv(results.rangesGSE, file="results.rangesGSE.csv")

# Perform gene set enrichment analysis of DMRs with missMethyl
gst.regionGSE <- goregion(results.rangesGSE, all.cpg=rownames(MvalGSE_no3), 
                          collection="GO", array.type="EPIC", plot.bias=TRUE,
                          sig.genes=F)
topGSA.gst.regionGSE<-topGSA(gst.regionGSE, n=Inf)
saveRDS(topGSA.gst.regionGSE,file="topGSA.gst.region.AFAB.GSE173717.4-5-23.RDS")

# Create a bubble plot of the top 10 GO terms for DMRs
topGSA.gst.regionGSE<-readRDS("topGSA.gst.region.AFAB.GSE173717.4-5-23.RDS")
top10_gst.regionGSE<-topGSA.gst.regionGSE %>% head(n=10)
ggplot(top10_gst.regionGSE, aes(x = DE, y = TERM, size = -log10(P.DE), 
                                label =(""))) +
    geom_point(alpha = 0.7, color = "blue") +
    geom_text(size = 5, vjust = 2) +
    scale_size(range = c(2, 10)) +
    labs(x = "Differentially Methylated Genes", y = "", 
         title = "Top 10 GO Terms for DMRs in Transmasculine Adults",
         size = "-log10(P-value)") +
    theme_bw() +
    theme(axis.text = element_text(size = 20), 
          axis.title.x=element_text(size = 20,
                                    margin = margin(t = 30, r = 0, b = 10, 
                                                    l = 0), 
                                    vjust = 2),
          axis.title.y=element_text(size = 20,
                                    margin = margin(t = 0, r = 152, b = 0, 
                                                    l = 10), 
                                    vjust = 2),
          legend.title=element_text(size = 16),
          legend.text=element_text(size = 14))+
    scale_x_continuous(limits = c(0, 19000))



#### Compare Direction of Methylation in Transmasculine Youth & Adults ####
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
# Load transmasculine youth GRanges from DMRcate
results.rangesAFAB<-readRDS("results.rangesAFAB.RDS")

setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/03-28-23")
# Load transmasculine adult GRanges from DMRcate
results.rangesGSE<-readRDS("results.rangesGSE.RDS")

setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

# Check which genes in youth DMRs are also in adult DMRs
youthDMR_in_adultDMR<-results.rangesAFAB[results.rangesAFAB$overlapping.genes 
                                         %in%results.rangesGSE$overlapping.genes,]
# ^ 3 of the 4 youth DMRs (20 genes) have genes that are also in adult DMRs
# ^ Only 1 DMR (1 gene) is statistically significant (with Fisher < 0.1)
write.csv(youthDMR_in_adultDMR, file="youthDMR_in_adultDMR.csv")

# Check which adult DMRs have the gene identified above
subset(results.rangesGSE, grepl("SLC1A4", overlapping.genes)) # in 2 DMRs
# ^ After 6 months, it is hypomethylated in both DMRs (meandiff > 0)

# Plot the youth DMR that overlaps SLC1A4
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")
groupAFAB<-readRDS("groupAFAB.RDS")
MvalAFAB_no3<-readRDS("MvalAFAB_no3.RDS")
cols <- c(1,2)[groupAFAB]
names(cols) <-groupAFAB
DMR.plot(ranges=results.rangesAFAB, dmr=2, CpGs=MvalAFAB_no3, phen.col=cols, 
         what="M", arraytype="EPIC", genome="hg19")
# ^ doesn't work; NAs?

# Plot the two adult DMRs (3479 and 56058) that overlap SLC1A4
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/GSE173717/03-28-23")
groupGSE<-readRDS("groupGSE.RDS")
MvalGSE_no3<-readRDS("MvalGSE_no3.RDS")
cols <- c(1,2)[groupGSE]
names(cols) <-groupGSE
DMR.plot(ranges=results.rangesGSE, dmr=3479, CpGs=MvalGSE_no3, phen.col=cols, 
         what="M", arraytype="EPIC", genome="hg19")
# ^ doesn't work; NAs?
DMR.plot(ranges=results.rangesGSE, dmr=56058, CpGs=MvalGSE_no3, phen.col=cols, 
         what="M", arraytype="EPIC", genome="hg19")
# ^ doesn't work; NAs?

# GRanges columns notes:
# seqnames = name of the chromosome where the interval's located
# ranges = start and end position of interval
# strand = forward (+), reverse (-), or strand-independent (*)
# no.cpgs = number of CpGs in DMR
# Fisher = multiple comparison statistic
# meandiff = difference bw means (Visit 1- Visit 2; negative number means
# Visit 2 is hypermethylated)



###### Compare DMR-Associated GO Terms in Transmasculine Youth & Adults ######

# Determine which youth DMRs overlap with GO terms containing "estrogen"
subset(topGSA.gst.regionAFAB, grepl("estrogen", TERM)) # 16; none significant
# Determine which adult DMRs overlap with GO terms containing "estrogen"
subset(topGSA.gst.regionGSE, grepl("estrogen", TERM)) # 16; none significant,
# ^ but two with 0.05 < P.DE < 0.1:
# ^ intracellular estrogen receptor signaling pathway P.DE = 0.064
# ^ nuclear estrogen receptor binding P.DE = 0.082

# Determine which youth DMRs overlap with GO terms containing "testosterone"
subset(topGSA.gst.regionAFAB, grepl("testosterone", TERM)) # 17; none significant
# Determine which adult DMRs overlap with GO terms containing "testosterone"
subset(topGSA.gst.regionGSE, grepl("testosterone", TERM)) # 17; none significant


# Select significant (p < 0.05) GO terms associated with youth DMRs
sigGO_DMR_AFAB<-topGSA.gst.regionAFAB %>% subset(., (P.DE < 0.05)) # all FDR = 1
# Select significant (p < 0.05) GO terms associated with adult DMRs
sigGO_DMR_GSE<-topGSA.gst.regionGSE %>% subset(., (P.DE < 0.05))

# Merge significant GO terms shared by youth and adult DMRs
siGO_DMR_youthadult<-merge(sigGO_DMR_AFAB,sigGO_DMR_GSE, by.x="TERM",
                           by.y="TERM")
# ^ P.DE.x = AFAB; P.DE.y = GSE

# Create a bubble plot of the GO terms significant for youth and adult DMRs
ggplot(siGO_DMR_youthadult) +
    geom_point(aes(x = DE.x, y = TERM, size = -log10(P.DE.x)),
               alpha = 0.7, color = "red") + # youth
    geom_point(aes(x = DE.y, y = TERM, size = -log10(P.DE.y)), 
               alpha = 0.7, color = "blue") + # adults
    labs(x = "Differentially Methylated Genes", y = "", 
         title = "Significant DMR GO Terms in Transmasculine Youth and Adults",
         size = "-log10(P-value)") +
    theme_bw()+
    theme(axis.text = element_text(size = 12), 
          axis.title.x=element_text(size = 20,
                                    margin = margin(t = 30, r = 0, b = 10, 
                                                    l = 0), 
                                    vjust = 2),
          axis.title.y=element_text(size = 20,
                                    margin = margin(t = 0, r = 240, b = 0,
                                                    l = 10), 
                                    vjust = 2),
          legend.title=element_text(size = 16),
          legend.text=element_text(size = 14))+
    scale_x_continuous(limits = c(0, 19000))





###### DMP and DMR Identification in Transmasculine Youth on Blockers ######
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

# Load objects from "pre-processing and normalization.R" code
# G954 metadata
G954sheet <- readRDS("G954sheet.RDS")
# G954 RGChannelSetExtended
rg <- readRDS("G954_EPIC_Data_0217rgsetext.RDS")
# G954 beta values
rg_betas <- readRDS("G954_EPIC_Data_0217rg_betas.RDS")

# Convert beta values to M-values using ENmix
MvalBlocker<-B2M(rg_betas)

# Select columns of interest (people only on blockers)
G954_Blocker<-G954sheet %>% subset(., (Sample_Name %in% 
                                           c("206839580019_R01C01",# 10002
                                             "206839580019_R02C01", # 10002_2
                                             "206839580032_R07C01",# 10017
                                             "206839580032_R08C01",# 10017_2
                                             "206839580066_R02C01",# 10018
                                             "206839580066_R01C01")))# 10018_2


# Create design matrix
groupBlocker <- factor(G954_Blocker$Visit,levels=c(1,2))
idBlocker <- factor(G954_Blocker$Subject.ID)
designBlocker <- model.matrix(~idBlocker + groupBlocker)
designBlocker

# Create character vector of Sample_Names of interest (people on blockers)
Blocker_string_use<-G954_Blocker %>% pull(Sample_Name)

# Select columns for samples of interest (people on blockers)
MvalBlocker_no3<-MvalBlocker[,(colnames(MvalBlocker) %in% Blocker_string_use)] 
saveRDS(MvalBlocker_no3,file="MvalBlocker_no3.RDS")

# Test for differential methylation using the lmFit and eBayes 
# functions from limma (input data = matrix of M-values)
fit.reducedBlocker <- lmFit(MvalBlocker_no3,designBlocker)
fit.reducedBlocker <- eBayes(fit.reducedBlocker, robust=TRUE)
# Check the numbers of hyper-methylated (1) and hypo-methylated (-1)
# using the decideTests function in limma 
summaryBlocker<-summary(decideTests(fit.reducedBlocker))
saveRDS(summaryBlocker,file="summaryBlocker.RDS")
decideTestsBlocker<-decideTests(fit.reducedBlocker)
saveRDS(decideTestsBlocker,file="decideTestsBlocker.RDS")
# Extract top 10 differentially methylated CpGs for Visit 1 versus Visit 2 
# using topTable
topBlocker<-topTable(fit.reducedBlocker,coef=4)
topBlocker
saveRDS(topBlocker,file="topBlocker.RDS")
# Set up the factor of interest
grpBlocker <- factor(G954_Blocker$Visit, labels=c(1,2))
# Extract Illumina negative control data
INCsBlocker <- getINCs(rg) 
# Select columns for samples of interest (people on blockers)
INCsBlocker_no3<-INCsBlocker[,(colnames(INCsBlocker) %in% Blocker_string_use)] 
saveRDS(INCsBlocker_no3,file="INCsBlocker_no3.RDS")
# Add negative control data to M-values
McBlocker <- rbind(MvalBlocker_no3,INCsBlocker_no3)
# Create vector marking negative controls in data matrix
ctl1Blocker <- rownames(McBlocker) %in% rownames(INCsBlocker_no3)
table(ctl1Blocker)
# Stage 1 analysis: perform RUV adjustment and fit to rank the CpGs with respect 
# to association with factor of interest
rfit1Blocker <- RUVfit(Y = McBlocker, X = grpBlocker, ctl = ctl1Blocker)
rfit2Blocker <- RUVadj(Y = McBlocker, fit = rfit1Blocker)
# Designate the CpGs that are LEAST associated with the factor of 
# interest (based on FDR-adjusted p-value) as empirical control probes (ECPs)
top1Blocker <- topRUV(rfit2Blocker, num=Inf, p.BH = 1)
ctl2Blocker <- rownames(MvalBlocker_no3) %in% 
    rownames(top1Blocker[top1Blocker$p.BH_X1.2 > 0.5,])
table(ctl2Blocker)
# Stage 2 analysis: use the ECPs to perform a second differential methylation 
# with RUV-inverse, which is adjusted for the unwanted variation estimated from 
# the data; perform RUV adjustment and fit
rfit3Blocker <- RUVfit(Y = MvalBlocker_no3, X = grpBlocker, ctl = ctl2Blocker)
# Error In function invvar: The fourth-smallest eigenvalue is more than 10 times 
#larger than the smallest eigenvalue.  This is not currently supported, and very 
#likely means something is wrong (perhaps a dimension was removed during preprocessing?).
# ^ this error can happen if number of variables >> sample size
rfit4Blocker <- RUVadj(Y = MvalBlocker_no3, fit = rfit3Blocker)
# Visualize top-ranked differentially methylated CpGs obtained from RUV 
# differential methylation analysis
    # p.BH = cutoff value for Benjamini-Hochberg adjusted p-values
    # p.BH_X1.2 = BH-adjusted p-values associated with factor of interest
fitdfBlocker<-as.data.frame(topRUV(rfit4Blocker,number=Inf))
saveRDS(fitdfBlocker,file="fitdfBlocker.RDS")
# Count p values < 0.05  
table(fitdfBlocker$p.BH_X1.2<0.05)

# Select all significant CpGs in dataset
top9357<-as.data.frame(topRUV(rfit4Blocker,number=9357))
saveRDS(top9357,file="top9357.RDS")
top9357<-readRDS("top9357.RDS")
sigCpGsBlocker<-rownames(top9357)
# Check number of genes to which significant CpGs are annotated
checkBlocker <- getMappedEntrezIDs(sig.cpg = sigCpGsBlocker)
length(checkBlocker$sig.eg)
# ^ ___ CpGs are differentially methylated in GAH pts between
# Visits 1 and 2, mapping to ___ different genes

# Perform gene set enrichment analysis of DMPs with missMethyl
gstBlocker <- gometh(sig.cpg=sigCpGsBlocker, all.cpg=rownames(fitdfBlocker), 
                     collection="GO",
                     plot.bias=TRUE, array.type="EPIC",sig.genes=T)
topGSA.gstBlocker<-topGSA(gstBlocker, n=Inf)
saveRDS(topGSA.gstBlocker,file="topGSA.gst.Blocker.sig.genes.RDS")


# Create a bubble plot of top 10 GO terms for DMPs
topGSA.gstBlocker<-readRDS("topGSA.gst.Blocker.sig.genes.RDS")
top10Blocker<-topGSA.gstBlocker %>% head(n=10)
ggplot(top10Blocker, aes(x = DE, y = TERM, size = -log10(P.DE), label =(""))) +
    geom_point(alpha = 0.7, color = "purple") +
    geom_text(size = 5, vjust = 2) +
    scale_size(range = c(2, 10)) +
    labs(x = "# Differentially Methylated Genes", y = "Gene Ontology Term",
         title = "Top 10 GO Terms in Transfeminine Youth on Blockers",
         size = "-log10(P-value)") +
    theme_bw()

# Annotate M-values with genomic information in DMRcate
myAnnotationBlocker <- cpg.annotate(object = MvalBlocker_no3, 
                                    datatype = "array", what = "M", 
                                    arraytype = c("EPIC"), 
                                    analysis.type = "differential", 
                                    design = designBlocker, 
                                    coef = 4)

# Identify DMRs; create DMResults object
DMRsBlocker <- dmrcate(myAnnotationBlocker, lambda=1000, C=2)
# Extract DMR range information; create GRanges object
results.rangesBlocker <- extractRanges(DMRsBlocker)
results.rangesBlocker
saveRDS(results.rangesBlocker, file="results.rangesBlocker.RDS")
# ^ ___ DMRs

# Perform gene set enrichment analysis of DMRs with missMethyl
gst.regionBlocker <- goregion(results.rangesBlocker, 
                              all.cpg=rownames(MvalBlocker_no3), 
                              collection="GO", array.type="EPIC", 
                              plot.bias=TRUE,
                              sig.genes=T
)
topGSA.gst.regionBlocker<-topGSA(gst.regionBlocker, n=Inf)
saveRDS(topGSA.gst.regionBlocker,file="topGSA.gst.region.Blocker.sig.genes.RDS")

# Create a bubble plot of the top 10 GO terms for DMRs
topGSA.gst.regionBlocker<-readRDS("topGSA.gst.region.Blocker.sig.genes.RDS")
top10.gst.regionBlocker<-topGSA.gst.regionBlocker %>% head(n=10)
ggplot(top10.gst.regionBlocker, aes(x = DE, y = TERM, size = -log10(P.DE), 
                                    label =(""))) +
    geom_point(alpha = 0.7, color = "purple") +
    geom_text(size = 5, vjust = 2) +
    scale_size(range = c(2, 10)) +
    labs(x = "# Differentially Methylated Genes", y = "Gene Ontology Term", 
         title = "Top 10 GO Terms for DMRs in Transfeminine Youth on Blockers",
         size = "-log10(P-value)") +
    theme_bw()


###### Transfeminine Youth on GAH ######
# Set working directory
setwd("/Volumes/nottermanlab/Gillian Hilscher/GA Study/Methylation Datasets/G954_EPIC_Data_0217/Processing/03-22-23")

# Load objects
# G954 metadata
G954sheet <- readRDS("G954sheet.RDS")
# G954 RGChannelSetExtended
rg <- readRDS("G954_EPIC_Data_0217rgsetext.RDS")
# G954 beta values
rg_betas <- readRDS("G954_EPIC_Data_0217rg_betas.RDS")

# Convert beta values to M-values using ENmix
MvalAMAB<-B2M(rg_betas)

# Select columns of interest (people transfeminine youth on GAH)
G954_AMAB<-G954sheet %>% subset(., (Sample_Name %in% 
                                        c("206842050024_R03C01", #10001
                                          "206842050024_R01C01", #10001_2
                                          "206839580032_R03C01",#10013
                                          "206839580032_R04C01", #10013_2
                                          "206839580066_R05C01", #10023
                                          "206839580066_R06C01"))) # 10023_2

# Create design matrix
groupAMAB <- factor(G954_AMAB$Visit,levels=c(1,2))
idAMAB <- factor(G954_AMAB$Subject.ID)
designAMAB <- model.matrix(~idAMAB + groupAMAB)
designAMAB

# Create character vector of Sample_Names of interest (transfeminine youth on 
# GAH)
AMAB_string_use<-G954_AMAB %>% pull(Sample_Name)

# Select columns for samples of interest (transfeminine youth on 
# GAH)
MvalAMAB_no3<-MvalAMAB[,(colnames(MvalAMAB) %in% AMAB_string_use)] 
saveRDS(MvalAMAB_no3,file="MvalAMAB_no3.RDS")

# Test for differential methylation using the lmFit and eBayes 
#functions from limma (input data = matrix of M-values)
fit.reducedAMAB <- lmFit(MvalAMAB_no3,designAMAB)
fit.reducedAMAB <- eBayes(fit.reducedAMAB, robust=TRUE)
# Check the numbers of hyper-methylated (1) and hypo-methylated (-1) 
# using the decideTests function in limma 
summaryAMAB<-summary(decideTests(fit.reducedAMAB))
saveRDS(summaryAMAB,file="summaryAMAB.RDS")
decideTestsAMAB<-decideTests(fit.reducedAMAB)
saveRDS(decideTestsAMAB,file="decideTestsAMAB.RDS")
# Extract the top 10 differentially methylated CpGs for Visit 1 versus Visit 2 
# using topTable
topAMAB<-topTable(fit.reducedAMAB,coef=4)
topAMAB
saveRDS(topAMAB,file="topAMAB.RDS")

# Set up the factor of interest
grpAMAB <- factor(G954_AMAB$Visit, labels=c(1,2))
# Extract Illumina negative control data
INCsAMAB <- getINCs(rg) 
# Select columns for samples of interest (transfeminine youth on GAH)
INCsAMAB_no3<-INCsAMAB[,(colnames(INCsAMAB) %in% AMAB_string_use)] 
saveRDS(INCsAMAB_no3,file="INCsAMAB_no3.RDS")
# Add negative control data to M-values
McAMAB <- rbind(MvalAMAB_no3,INCsAMAB_no3)
# Create vector marking negative controls in data matrix
ctl1AMAB <- rownames(McAMAB) %in% rownames(INCsAMAB_no3)
table(ctl1AMAB)
# Stage 1 analysis: perform RUV adjustment and fit to rank the CpGs with respect 
# to association with factor of interest
rfit1AMAB <- RUVfit(Y = McAMAB, X = grpAMAB, ctl = ctl1AMAB) # Stage 1 analysis
rfit2AMAB <- RUVadj(Y = McAMAB, fit = rfit1AMAB)
# Designate the CpGs that are LEAST associated with the factor of 
# interest (based on FDR-adjusted p-value) as empirical control probes (ECPs)
top1AMAB <- topRUV(rfit2AMAB, num=Inf, p.BH = 1)
ctl2AMAB <- rownames(MvalAMAB_no3) %in% 
    rownames(top1AMAB[top1AMAB$p.BH_X1.2 > 0.5,])
table(ctl2AMAB)
# Stage 2 analysis: use the ECPs to perform a second differential methylation 
# with RUV-inverse, which is adjusted for the unwanted variation estimated from 
# the data; perform RUV adjustment and fit
rfit3AMAB <- RUVfit(Y = MvalAMAB_no3, X = grpAMAB, ctl = ctl2AMAB)
rfit4AMAB <- RUVadj(Y = MvalAMAB_no3, fit = rfit3AMAB)
# Visualize top-ranked differentially methylated CpGs obtained from RUV 
# differential methylation analysis
    # p.BH = cutoff value for Benjamini-Hochberg adjusted p-values
    # p.BH_X1.2 = BH-adjusted p-values associated with factor of interest
fitdfAMAB<-as.data.frame(topRUV(rfit4AMAB,number=Inf))
saveRDS(fitdfAMAB,file="fitdfAMAB.RDS")
# Count p values < 0.05 
table(fitdfAMAB$p.BH_X1.2<0.05) # just look at p_X1.2 instead?
# ^ none of the BH adjusted p values are less than 0.05 (no significant DMPs)

# Select all significant CpGs in dataset
top9357<-as.data.frame(topRUV(rfit4AMAB,number=9357))
saveRDS(top9357,file="top9357.RDS")
top9357<-readRDS("top9357.RDS")
sigCpGsAMAB<-rownames(top9357)
# Check number of genes to which significant CpGs are annotated
checkAMAB <- getMappedEntrezIDs(sig.cpg = sigCpGsAMAB)
length(checkAMAB$sig.eg)
# ^ ___ CpGs are differentially methylated in GAH pts between
# Visits 1 and 2, mapping to ___ different genes

# Perform gene set enrichment analysis of DMPs with missMethyl
gstAMAB <- gometh(sig.cpg=sigCpGsAMAB, all.cpg=rownames(fitdfAMAB), 
                  collection="GO",
                  plot.bias=TRUE, array.type="EPIC",sig.genes=T)
topGSA.gstAMAB<-topGSA(gstAMAB, n=Inf)
saveRDS(topGSA.gstAMAB,file="topGSA.gst.AMAB.sig.genes.RDS")


# Create a bubble plot of top 10 GO terms for DMPs
topGSA.gstAMAB<-readRDS("topGSA.gst.AMAB.sig.genes.RDS")
top10AMAB<-topGSA.gstAMAB %>% head(n=10)
ggplot(top10AMAB, aes(x = DE, y = TERM, size = -log10(P.DE), label =(""))) +
    geom_point(alpha = 0.7, color = "orange") +
    geom_text(size = 5, vjust = 2) +
    scale_size(range = c(2, 10)) +
    labs(x = "# Differentially Methylated Genes", y = "Gene Ontology Term", 
         title = "Top 10 GO Terms in Transfeminine Youth on GAH",
         size = "-log10(P-value)") +
    theme_bw()

# Annotate M-values with genomic information in DMRcate
myAnnotationAMAB <- cpg.annotate(object = MvalAMAB_no3, datatype = "array", 
                                 what = "M", 
                                 arraytype = c("EPIC"), 
                                 analysis.type = "differential", 
                                 design = designAMAB, 
                                 coef = 4)
# Identify DMRs; create DMResults object
DMRsAMAB <- dmrcate(myAnnotationAMAB, lambda=1000, C=2)
# ^ Error: The FDR you specified in cpg.annotate() returned no significant CpGs, 
# hence there are no DMRs.  Try specifying a value of 'pcutoff' in dmrcate() 
# and/or increasing 'fdr' in cpg.annotate().

# Extract DMR range information; create GRanges object
results.rangesAMAB <- extractRanges(DMRsAMAB)
results.rangesAMAB
saveRDS(results.rangesAMAB, file="results.rangesAMAB.RDS")
# ^ ___ DMRs

# Perform gene set enrichment analysis of DMRs with missMethyl
gst.regionAMAB <- goregion(results.rangesAMAB, all.cpg=rownames(MvalAMAB_no3), 
                           collection="GO", array.type="EPIC", plot.bias=TRUE,
                           sig.genes=T
)
topGSA.gst.regionAMAB<-topGSA(gst.regionAMAB, n=Inf)
saveRDS(topGSA.gst.regionAMAB,file="topGSA.gst.region.AMAB.sig.genes.RDS")

# Create a bubble plot of the top 10 GO terms for DMRs
topGSA.gst.regionAMAB<-readRDS("topGSA.gst.region.AMAB.sig.genes.RDS")
top10.gst.regionAMAB<-topGSA.gst.regionAMAB %>% head(n=10)
ggplot(top10.gst.regionAMAB, aes(x = DE, y = TERM, size = -log10(P.DE), 
                                 label =(""))) +
    geom_point(alpha = 0.7, color = "purple") +
    geom_text(size = 5, vjust = 2) +
    scale_size(range = c(2, 10)) +
    labs(x = "# Differentially Methylated Genes", y = "Gene Ontology Term", 
         title = "Top 10 GO Terms for DMRs in Transfeminine Youth on GAH",
         size = "-log10(P-value)") +
    theme_bw()