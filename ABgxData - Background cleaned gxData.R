require("ggplot2")
require("tidyverse")
require("dplyr")
require("tidyr")
require("pcaMethods")
require("Biobase")
require("limma")
require("biomaRt")
require("gridExtra")
require("viridis")
library(WGCNA)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
Ensdb = EnsDb.Hsapiens.v86
require("Rgraphviz")
require("DOSE")
require("ReactomePA")
require("pathview")
require("GOSemSim")
require("enrichplot")
require("clusterProfiler")
require("pathfindR")
require("RDAVIDWebService")
require("GOplot")
require("topGO")
require("goseq")

options(stringsAsFactors = FALSE)

#-----------------------------------------------------------------------------#
#1: Importing the data and inspecting sample information
#-----------------------------------------------------------------------------

#1a - Importing data

setwd("D:/Documents/Github/Research_project_1/CardiomyopathyGroup2") # setting working directory
gxData <- read.delim("D:/Documents/Github/Research_project_1/MAGNET_GeneExpressionData_CPM_19112020.txt", row.names = 1, header = TRUE) # importing gene expression data
sampleInfo <- read.delim("D:/Documents/Github/Research_project_1/MAGNET_SampleData_19112020.txt", row.names = 1, header = TRUE) # importing phenotype data
geneTotExonLengths <- read.delim("D:/Documents/Github/Research_project_1/MAGNET_exonLengths.txt", as.is = T, row.names = 1) # importing gene length data

# Transforming gxData to FPKM from CPM

all(rownames(geneTotExonLengths) == rownames(gxData))
cpm2fpkm <- function(x) {
  geneTotExonLengths_kb <- geneTotExonLengths[, 1] / 1E3
  .t <- 2^(x) / geneTotExonLengths_kb
}
gxData_fpkm <- cpm2fpkm(gxData)

#setting ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#-------------------------------------------------------------------------------
# 1 - Annotating genes with entrez IDs and biotype.  
#-------------------------------------------------------------------------------
EntrezIDs = bitr(row.names(gxData), fromType = "ENSEMBL",
                 toType = c("ENTREZID"),
                 OrgDb = org.Hs.eg.db)
EntrezIDs = distinct(EntrezIDs, ENSEMBL, .keep_all=TRUE)

#n = data.frame(table(EntrezIDs$ENSEMBL))
#n[n$Freq!=1,]
#dim(EntrezIDs[n$Freq>1,])


#subsetting all genes that have Entrez IDs
gxEntrez = gxData[EntrezIDs[,1],]
geneInfo = getBM(attributes = c("gene_biotype", "ensembl_gene_id"), 
                 filters = c("ensembl_gene_id"), 
                 values = row.names(gxEntrez), 
                 mart = ensembl
)
# checking for missing/repeated values. - there are none. 
n = data.frame(table(geneInfo$ensembl_gene_id))
geneInfo[n$Freq!=1,]

#-------------------------------------------------------------------------------
# 3. Checking background noise level. 
#-------------------------------------------------------------------------------

# Finding all Y-chromosome genes
yGeneID <- getBM(
  attributes = c("ensembl_gene_id"), 
  filter = c("chromosome_name", "ensembl_gene_id"),
  values = list("Y", row.names(gxData_fpkm)), mart = ensembl)

# Subsetting gene expression data for females and only genes on Y chromosome.

yGene <- gxData_fpkm[yGeneID[, 1], sampleInfo$Sex == "Female"]

# Finding the Mean expression of these genes. Since CPM has been normalized for 
# between-sample comparisons, and FPKM is derived from this CPM, it is possible
# to use FPKM for between-sample comparisons as well. 

background <- mean(rowMeans(yGene)) 
rmeans <- as.matrix(rowMeans(gxData_fpkm)) # Calculating average expression across all individuals for all genes.
colnames(rmeans)<-c('Row Means')

# Subsetting genes that are expressed above background level 

aboveBackground<- NULL # Creating empty variable
belowBackground<- NULL
for (i in 1:nrow(gxData_fpkm)) {
  if ((rmeans[i]) > background) {
    aboveBackground <- rbind(aboveBackground, (row.names(gxData_fpkm)[i])) # adding gene ID to the variable
  }
  if ((rmeans[i])<=background){
    belowBackground <- rbind(belowBackground, (row.names(gxData_fpkm)[i]))
  }
}
colnames(aboveBackground)<-c("Gene ID")
colnames(belowBackground)<-c("Gene ID")

#Cleaned, above background data-  to be used for further analysis. 
ABgxData = gxData[aboveBackground,]

ABgeneTotExonLengths = geneTotExonLengths[aboveBackground,]

write.table(ABgxData, file = "D:/Documents/Github/Research_project_1/ABgxData.txt", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(ABgeneTotExonLengths, file = "D:/Documents/Github/Research_project_1/ABgeneTotExonLengths.txt", sep = "\t", row.names = TRUE, quote = FALSE)

#---------------------------------------------------------------------------------
# 3 - Differential gene expression analysis
#---------------------------------------------------------------------------------
# Creating an expressionSet object
eset <- ExpressionSet(as.matrix(ABgxData),
                      phenoData = AnnotatedDataFrame(ABsampleInfo),
                      featureData = AnnotatedDataFrame(ABgeneTotExonLengths)
)

# Creating a design matrix with no intercept and 'Disease', 'Sex' and 'Ethnicity' 
# as explanatory variables. All covariates were included in the linear model to 
# account for any clustering that may be present in the dataset which may be 
# visible only in 5th or later principle components. However, only disease
# was seen to the reason for clustering in the first 5 principle components. 

#con = list(relevel(factor(sampleInfo$Disease), ref = "Donor"))
#names(con[[1]])=rownames(sampleInfo)

design <- model.matrix(~ 0 + Disease + Sex + Ethnicity ,data = pData(eset))


# Defining contrasts
cont2 <- makeContrasts(
  DCMvDon = DiseaseDCM - DiseaseDonor,
  HCMvDon = DiseaseHCM - DiseaseDonor,
  PPCMvDon = DiseasePPCM - DiseaseDonor, 
  DCMvHCM = DiseaseDCM - DiseaseHCM, 
  DCMvPPCM = DiseaseDCM - DiseasePPCM, 
  HCMvPPCM = DiseaseHCM - DiseasePPCM, levels = design
)


# Fitting linear model
fit <- lmFit(eset, design)

# Fitting the defined contrasts
fit2 <- contrasts.fit(fit, cont2)

# eBayes analysis
fit2 <- eBayes(fit2, trend=T)

# Classifying genes as up, down or not significantly regulated.
results <- decideTests(fit2)
summary(results)
t4<-summary(results)

# Differentially expressed genes in DCM vs Donors- logFC threshold was kept to log2(1.5).
# topTable and topTreat both use lfc in log2. So setting lfc to log2(1.5) would mean lfc~0.585
# which means expression has to be at least 1.5 times to be over/underexpressed.
# Benjamini-Hochberg correction was applied to p value. 0.05 was set as the cutoff
# for this adjusted p value (False Discovery Rate).

degDCM <- topTable(fit2, coef = "DCMvDon", adjust.method = "BH", lfc = log2(1.5),
                   number = nrow(fit2), sort.by = "none", p.value = 0.05
)
degHCM <- topTable(fit2, coef = "HCMvDon", adjust.method = "BH", lfc = log2(1.5),
                   number = nrow(fit2), sort.by = "none", p.value = 0.05
)
degPPCM <- topTable(fit2, coef = "PPCMvDon", adjust.method = "BH", lfc = log2(1.5),
                    number = nrow(fit2), sort.by = "none", p.value = 0.05
)
degDCMvHCM <- topTable(fit2, coef = "DCMvHCM", adjust.method = "BH", lfc = log2(1.5),
                       number = nrow(fit2), sort.by = "none", p.value = 0.05
)

# Subsetting up and down regulated, and all differentially expressed genes in DCM individuals.
upRegulatedGenes<-degDCM[degDCM$logFC>0,] 
downRegulatedGenes<-degDCM[degDCM$logFC<0,]

# Histogram of p values as a sanity check
hist(degDCM[, "adj.P.Val"], main="Adjusted P Value")

# Finding genes that are DE in all three disease groups
deg_overlap <- row.names(ABgxData[row.names(ABgxData) %in% row.names(degDCM) & row.names(ABgxData) %in% row.names(degHCM) & row.names(ABgxData) %in% row.names(degPPCM),])
overlapDEgenes = ABgxData[deg_overlap,]
ppcmDEgenes = ABgxData[row.names(degPPCM),]
dcmDEgenes = ABgxData[row.names(degDCM),]
hcmDEgenes = ABgxData[row.names(degHCM),]

