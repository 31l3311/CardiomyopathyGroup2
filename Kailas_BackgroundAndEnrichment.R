require("dplyr")
require("Biobase")
require("limma")
require("biomaRt")
require("xlsx")
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
Ensdb <- EnsDb.Hsapiens.v86
require("enrichplot")
require("clusterProfiler")
require("pathfindR")
require("topGO")

options(stringsAsFactors = FALSE)

#-----------------------------------------------------------------------------#
# 1: Importing the data and inspecting sample information
#-----------------------------------------------------------------------------
# Importing data

setwd("C:/Users/kaila/Documents/MSc/Period 2/Experimental methods and data management/R skill sessions/Data") # setting working directory
gxData <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", row.names = 1, header = TRUE) # importing gene expression data
sampleInfo <- read.delim("MAGNET_SampleData_19112020.txt", row.names = 1, header = TRUE) # importing phenotype data
geneTotExonLengths <- read.delim("MAGNET_exonLengths.txt", as.is = T, row.names = 1) # importing gene length data

# Transforming gxData to FPKM from CPM

all(rownames(geneTotExonLengths) == rownames(gxData))
cpm2fpkm <- function(x) {
  geneTotExonLengths_kb <- geneTotExonLengths[, 1] / 1E3
  .t <- 2^(x) / geneTotExonLengths_kb
}
gxData_fpkm <- cpm2fpkm(gxData)

# setting ensembl database
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#-------------------------------------------------------------------------------
# 2. Checking background noise level.
#-------------------------------------------------------------------------------

YGenes <- getBM(
  attributes = "ensembl_gene_id",
  filters = c("chromosome_name", "ensembl_gene_id"),
  values = list("Y", row.names(gxData)),
  mart = ensembl
)
YGenes <- unlist(c(YGenes))

# Calculating median noise gene
i <- 0
medians <- medianGenes <- c()
for (s in colnames(gxData_fpkm[, which(sampleInfo$Sex == "Female")])) {
  i <- i + 1
  medians[i] <- median(gxData_fpkm[YGenes, s])
  medianGenes[i] <- YGenes[which.min(abs(gxData_fpkm[YGenes, s] - medians[i]))] # Here is the line calculating the Y-gene with median expression per female sample
}
names(medians) <- colnames(gxData_fpkm[, which(sampleInfo$Sex == "Female")])
names(medianGenes) <- colnames(gxData_fpkm[, which(sampleInfo$Sex == "Female")])
medianNoiseGene <- names(table(medianGenes))[as.vector(table(medianGenes)) == max(table(medianGenes))]

threshold <- apply(gxData[medianNoiseGene, which(sampleInfo$Sex == "Female")], 1, FUN = mean)

# Determining genes expressed above background based on CEFIC rules.
# Threshold determined above (~5.5 cpm). - A gene is considered lowly expressed
# and left out of analysis if more than 25 percent of samples in each of the 
# disease groups have an expression less than threshold. 
# This correction also accounts for 'spikes', i.e. one large expression value, 
# highly influences the average. 

matCounts <- as.matrix(gxData)
diseases <- c("DCM", "HCM", "PPCM")
qualityTest <- NULL
for (gene in row.names(gxData)) {
  Check3 <- c()
  for (s in diseases) {
    Samples <- colnames(gxData[, which(sampleInfo[, 3] == s)])
    matCounts2 <- matCounts[gene, Samples]
    Check1 <- sum(matCounts2 >= threshold) >= 0.75 * length(Samples)   
    one <- max(matCounts2) - median(matCounts2)                        #checking for spikes
    two <- sum(matCounts2) / (length(Samples) - 1)
    Check2 <- one < two
    Check3 <- c(Check3, Check1 & Check2)
  }
  qualityTest <- c(qualityTest, sum(Check3)>= 1)
}

# Subsetting genes that are expressed above background level
ABgxData <- gxData[qualityTest, ]
ABgeneTotExonLengths <- data.frame(geneTotExonLengths[qualityTest, ])
row.names(ABgeneTotExonLengths) <- row.names(ABgxData)
colnames(ABgeneTotExonLengths)[1] <- "ExonLegthBp"

#-------------------------------------------------------------------------------
# 3 - Annotating genes with entrez IDs and biotype
#-------------------------------------------------------------------------------
EntrezIDs <- bitr(row.names(ABgxData),
                  fromType = "ENSEMBL", # this does not give a 1-to-1 mapping. Many repeats.
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db
)
EntrezIDs <- distinct(EntrezIDs, ENSEMBL, .keep_all = TRUE) # Keep only the first instance of each ID
row.names(EntrezIDs) <- EntrezIDs[, 1]

# subsetting all genes that have Entrez IDs
ABgxEntrez <- ABgxData[EntrezIDs[, 1], ]

#biotype <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
#                 filters = "ensembl_gene_id", 
#                 values = row.names(ABgxData),
#                 mart = ensembl)
write.table(ABgxData, file = "ABgxData.txt", sep="\t")
save(sampleInfo, YGenes, ABgxData, ABgeneTotExonLengths, EntrezIDs, file = "DE analysis - Input Data.RData")
#---------------------------------------------------------------------------------
# 4 - Differential gene expression analysis
#---------------------------------------------------------------------------------
# Load saved data 
lnames <- load(file = "DE analysis - Input Data.RData")
# The variable lnames contains the names of loaded variables
lnames

# Creating an expressionSet object
eset <- ExpressionSet(as.matrix(ABgxData),
                      phenoData = AnnotatedDataFrame(sampleInfo),
                      featureData = AnnotatedDataFrame(ABgeneTotExonLengths)
)

# Creating a design matrix with no intercept and 'Disease', 'Sex' and 'Ethnicity'
# as explanatory variables. All covariates were included in the linear model to
# account for any clustering that may be present in the dataset which may be
# visible only in 5th or later principle components. However, only disease
# was seen to the reason for clustering in the first 5 principle components.

# Creating a new variable for contrasting sex and disease
disease.sex <- paste(sampleInfo$Disease, sampleInfo$Sex, sep = ".")

# Creating design matrix
design <- model.matrix(~ 0 + disease.sex + Ethnicity, data = pData(eset))

# Creating contrasts
cont2 <- makeContrasts(
  DCMvDon = (disease.sexDCM.Female + disease.sexDCM.Male) - (disease.sexDonor.Female + disease.sexDonor.Male),
  HCMvDon = (disease.sexHCM.Female + disease.sexHCM.Male) - (disease.sexDonor.Female + disease.sexDonor.Male),
  PPCMvDon = disease.sexPPCM.Female - disease.sexDonor.Female,
  InteractionDCM = (disease.sexDCM.Female - disease.sexDCM.Male) - (disease.sexDonor.Female - disease.sexDonor.Male),
  InteractionHCM = (disease.sexHCM.Female - disease.sexHCM.Male) - (disease.sexDonor.Female - disease.sexDonor.Male),
  levels = design
)

# Fitting linear model
fit <- lmFit(eset, design)

# Fitting the defined contrasts
fit2 <- contrasts.fit(fit, cont2)

# eBayes analysis
fit2 <- eBayes(fit2, trend = T)

# Classifying genes as up, down or not significantly regulated.
results <- decideTests(fit2)
summary(results)

# Significantly differentially expressed genes - logFC threshold was kept to log2(1.5).
# topTable uses log2 fold change. Setting lfc to log2(1.5) would mean lfc~0.585
# which means expression has to be at least 1.5 times to be over/underexpressed.
# Benjamini-Hochberg correction was applied to p value. 0.05 was set as the cutoff
# for this adjusted p value (False Discovery Rate).

resDCM <- topTable(fit2,
                   coef = "DCMvDon", adjust.method = "BH", lfc = log2(1.5),
                   number = nrow(fit2), sort.by = "logFC", p.value = 0.05
)
resHCM <- topTable(fit2,
                   coef = "HCMvDon", adjust.method = "BH", lfc = log2(1.5),
                   number = nrow(fit2), sort.by = "logFC", p.value = 0.05
)
resPPCM <- topTable(fit2,
                    coef = "PPCMvDon", adjust.method = "BH", lfc = log2(1.5),
                    number = nrow(fit2), sort.by = "logFC", p.value = 0.05
)
resIntDCM <- topTable(fit2,
                      coef = "InteractionDCM", adjust.method = "BH", lfc = log2(1.5),
                      number = nrow(fit2), sort.by = "logFC", p.value = 0.05
)
resIntHCM <- topTable(fit2,
                      coef = "InteractionHCM", adjust.method = "BH", lfc = log2(1.5),
                      number = nrow(fit2), sort.by = "logFC", p.value = 0.05
)
# Subsetting up and down regulated, and all differentially expressed genes in diseased individuals.

#upRegulatedDCM <- resDCM[resDCM$logFC > 0, ]
#downRegulatedDCM <- resDCM[resDCM$logFC < 0, ]

#upRegulatedHCM <- resHCM[resHCM$logFC > 0, ]
#downRegulatedHCM <- resHCM[resHCM$logFC < 0, ]

#upRegulatedPPCM <- resPPCM[resDCM$logFC > 0, ]
#downRegulatedPPCM <- resPPCM[resPPCM$logFC < 0, ]

# Histogram of p values
#hist(resDCM[, "adj.P.Val"], main = "Adjusted P Value for DCM")
#hist(resHCM[, "adj.P.Val"], main = "Adjusted P Value for HCM")
#hist(resPPCM[, "adj.P.Val"], main = "Adjusted P Value for PPCM")
#hist(resIntDCM[, "adj.P.Val"], main = "Adjusted P Value for Interaction DCM")
#hist(resIntHCM[, "adj.P.Val"], main = "Adjusted P Value for Interaction HCM")

# listing genes that are DE in all three disease groups
deg_overlap <- row.names(ABgxData[row.names(ABgxData) %in% row.names(resDCM) & row.names(ABgxData) %in% row.names(resHCM) & row.names(ABgxData) %in% row.names(resPPCM), ])

# subsetting DE genes by disease
#overlapDEgenes <- ABgxData[deg_overlap, ]
# ppcmDEgenes = ABgxData[row.names(resPPCM),]
# dcmDEgenes = ABgxData[row.names(resDCM),]
# hcmDEgenes = ABgxData[row.names(resHCM),]

# Listing all significantly differentially expressed genes. 
# degList <- row.names(ABgxData) %in% row.names(resDCM) | row.names(ABgxData) %in% row.names(resHCM) | row.names(ABgxData) %in% row.names(resPPCM)
# ABgxDataDEG <- ABgxData[degList, ]

#-------------------------------------------------------------------------------
# Creating dataframe of differentially expressed genes for WGCNA modules
#-------------------------------------------------------------------------------
# Function to make a dataframe from list with vectors of different lengths
listToDF <- function(list){
  sapply(list, "length<-", max(lengths(list)))
}

DEgenes <- list(row.names(resDCM), row.names(resHCM), row.names(resPPCM), row.names(resIntDCM), row.names(resIntHCM))
names(DEgenes) <- c("DEgenes_DCM", "DEgenes_HCM", "DEgenes_PPCM", "DEgenes_IntDCM","DEgenes_IntHCM")

# Create a data frame with DEGs, one column per disease
DEgenes<-listToDF(DEgenes)

write.table(DEgenes, file="allDEgenes.txt", sep='\t', row.names = FALSE)

#-------------------------------------------------------------------------------
# 6. Preparing data for ORA/pathway analysis
#
# Creating dataframe for each contrast with Ensembl/Entrez gene id, logFC, adj.p.value
# USING SIGNIFICANT DEGs only.
#-------------------------------------------------------------------------------
# Ensembl
# significant DEGs for dcm
DEGdcmEnsembl <- resDCM[, c(2, 6)]
DEGdcmEnsembl[, 3] <- row.names(DEGdcmEnsembl)
DEGdcmEnsembl <- DEGdcmEnsembl[, c(3, 1, 2)]
colnames(DEGdcmEnsembl)[1] <- "EnsemblID"

# significant DEGs for hcm
DEGhcmEnsembl <- resHCM[, c(2, 6)]
DEGhcmEnsembl[, 3] <- row.names(DEGhcmEnsembl)
DEGhcmEnsembl <- DEGhcmEnsembl[, c(3, 1, 2)]
colnames(DEGhcmEnsembl)[1] <- "EnsemblID"

# significant DEGs for ppcm
DEGppcmEnsembl <- resPPCM[, c(2, 6)]
DEGppcmEnsembl[, 3] <- row.names(DEGppcmEnsembl)
DEGppcmEnsembl <- DEGppcmEnsembl[, c(3, 1, 2)]
colnames(DEGppcmEnsembl)[1] <- "EnsemblID"

# significant DEGs for Interaction DCM
DEGdcmINTEnsembl <- resIntDCM[, c(2, 6)]
DEGdcmINTEnsembl[, 3] <- row.names(DEGdcmINTEnsembl)
DEGdcmINTEnsembl <- DEGdcmINTEnsembl[, c(3, 1, 2)]
colnames(DEGdcmINTEnsembl)[1] <- "EnsemblID"

# significant DEGs for Interaction HCM
DEGhcmINTEnsembl <- resIntHCM[, c(2, 6)]
DEGhcmINTEnsembl[, 3] <- row.names(DEGhcmINTEnsembl)
DEGhcmINTEnsembl <- DEGhcmINTEnsembl[, c(3, 1, 2)]
colnames(DEGhcmINTEnsembl)[1] <- "EnsemblID"

# Entrez
# dataframe of significant DEGs DCM (entrez id, logfc, adjusted pvalue)
DEGdcmEntrez <- DEGdcmEnsembl
DEGdcmEntrez[, 1] <- EntrezIDs[DEGdcmEntrez[, 1], 2]
colnames(DEGdcmEntrez)[1] <- "ENTREZID"
DEGdcmEntrez <- DEGdcmEntrez[!is.na(DEGdcmEntrez[, 1]), ]
DEGdcmEntrez <- DEGdcmEntrez[order(DEGdcmEntrez$logFC, decreasing = TRUE), ]

# dataframe of significant DEGs HCM (entrez id, logfc, adjusted pvalue)
DEGhcmEntrez <- DEGhcmEnsembl
DEGhcmEntrez[, 1] <- EntrezIDs[DEGhcmEntrez[, 1], 2]
colnames(DEGhcmEntrez)[1] <- "ENTREZID"
DEGhcmEntrez <- DEGhcmEntrez[!is.na(DEGhcmEntrez[, 1]), ]
DEGhcmEntrez <- DEGhcmEntrez[order(DEGhcmEntrez$logFC, decreasing = TRUE), ]

# dataframe of significant DEGs PPCM (entrez id, logfc, adjusted pvalue)
DEGppcmEntrez <- DEGppcmEnsembl
DEGppcmEntrez[, 1] <- EntrezIDs[DEGppcmEntrez[, 1], 2]
colnames(DEGppcmEntrez)[1] <- "ENTREZID"
DEGppcmEntrez <- DEGppcmEntrez[!is.na(DEGppcmEntrez[, 1]), ]
DEGppcmEntrez <- DEGppcmEntrez[order(DEGppcmEntrez$logFC, decreasing = TRUE), ]

# dataframe of significant DEGs InteractionDCM (entrez id, logfc, adjusted pvalue)
DEGdcmINTEntrez <- DEGdcmINTEnsembl
DEGdcmINTEntrez[, 1] <- EntrezIDs[DEGdcmINTEntrez[, 1], 2]
colnames(DEGdcmINTEntrez)[1] <- "ENTREZID"
DEGdcmINTEntrez <- DEGdcmINTEntrez[!is.na(DEGdcmINTEntrez[, 1]), ]
DEGdcmINTEntrez <- DEGdcmINTEntrez[order(DEGdcmINTEntrez$logFC, decreasing = TRUE), ]

# dataframe of significant DEGs InteractionHCM (entrez id, logfc, adjusted pvalue)
DEGhcmINTEntrez <- DEGhcmINTEnsembl
DEGhcmINTEntrez[, 1] <- EntrezIDs[DEGhcmINTEntrez[, 1], 2]
colnames(DEGhcmINTEntrez)[1] <- "ENTREZID"
DEGhcmINTEntrez <- DEGhcmINTEntrez[!is.na(DEGhcmINTEntrez[, 1]), ]
DEGhcmINTEntrez <- DEGhcmINTEntrez[order(DEGhcmINTEntrez$logFC, decreasing = TRUE), ]
#-------------------------------------------------------------------------------
# Creating dataframe for each disease of all gene IDs(ensembl), logfc, adj.p.value
#-------------------------------------------------------------------------------
adj.p <- list()
foldChange <- list()
geneName <- list()
for (i in 1:5) {
  x <- topTable(fit2, coef = i, adjust.method = "BH", sort.by = "logFC", number = nrow(fit2), p.value = 1)
  foldChange[[i]] <- x$logFC
  adj.p[[i]] <- x$adj.P.Val
  geneName[[i]] <- row.names(x)
}
foldChange <- data.frame(
  geneName[[1]], foldChange[[1]], adj.p[[1]],
  geneName[[2]], foldChange[[2]], adj.p[[2]],
  geneName[[3]], foldChange[[3]], adj.p[[3]],
  geneName[[4]], foldChange[[4]], adj.p[[4]],
  geneName[[5]], foldChange[[5]], adj.p[[5]]
)
colnames(foldChange) <- c(
  "EnsemblID", "DCMfc", "adjDCMp",
  "EnsemblID", "HCMfc", "adjHCMp",
  "EnsemblID", "PPCMfc", "adjPPCMp",
  "EnsemblID", "INTdcmFC", "adjINTDCMp",
  "EnsemblID", "INThcmFC", "adjINTHCMp"
)
rm(x)

allgenes_DCM <- foldChange[, c(1, 2, 3)]
allgenes_HCM <- foldChange[, c(4, 5, 6)]
allgenes_PPCM <- foldChange[, c(7, 8, 9)]
allgenes_INTdcm <- foldChange[, c(10, 11, 12)]
allgenes_INThcm <- foldChange[, c(13, 14, 15)]

# Creating .txt files of allgenes_ dataframes for any online GSEA tools that may be used
write.table(allgenes_DCM,file = "allgenes_DCM.csv", sep=",")
write.table(allgenes_HCM,file = "allgenes_HCM.txt", sep="\t")
write.table(allgenes_PPCM,file = "allgenes_PPCM.txt", sep="\t")
write.table(allgenes_INTdcm,file = "allgenes_INTdcm.txt", sep="\t")
write.table(allgenes_INThcm,file = "allgenes_INThcm.txt", sep="\t")
#-------------------------------------------------------------------------------
# Annotating the allgenes_ dataframes with gene symbols
#-------------------------------------------------------------------------------
diseases <- c("DCM", "HCM", "PPCM", "INTdcm", "INThcm")
for (i in diseases) {
  source <- get(paste("allgenes_", as.character(i), sep = ""))
  createdVar <- paste("inputList", as.character(i), sep = "_")
  assign(
    createdVar,
    bitr(source[, 1], fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
  )
}
#-------------------------------------------------------------------------------
# Creating inputLists for pathfindR - input must be a data frame with columns 
# containing: Gene Symbols, Change Values (optional) and p values.
#-------------------------------------------------------------------------------
# DCM
row.names(inputList_DCM) <- NULL
# inputList_DCM = inputList_DCM[order(inputList_DCM$ENSEMBL, decreasing = F),]
# inputList_DCM = inputList_DCM[!is.na(inputList_DCM[,1])&!is.na(inputList_DCM[,2]),]
# inputList_DCM = inputList_DCM[-c(443,909),]#after manual checking, removing wrong annotation
inputList_DCM <- distinct(inputList_DCM, ENSEMBL, .keep_all = TRUE)
inputList_DCM <- distinct(inputList_DCM, SYMBOL, .keep_all = TRUE)
inputList_DCM[, c(3, 4)] <- (allgenes_DCM[allgenes_DCM$EnsemblID %in% inputList_DCM[, 1], c(2, 3)])
inputList_DCM <- data.frame(inputList_DCM[, c(2, 3, 4)])

# HCM
row.names(inputList_HCM) <- NULL
# inputList_HCM = inputList_HCM[!is.na(inputList_HCM[,1])&!is.na(inputList_HCM[,2]),]
# inputList_HCM = inputList_HCM[-c(1419),]#after manual checking, removing wrong annotation
inputList_HCM <- distinct(inputList_HCM, ENSEMBL, .keep_all = TRUE)
inputList_HCM <- distinct(inputList_HCM, SYMBOL, .keep_all = TRUE)
inputList_HCM[, c(3, 4)] <- (allgenes_HCM[allgenes_HCM$EnsemblID %in% inputList_HCM[, 1], c(2, 3)])
inputList_HCM <- data.frame(inputList_HCM[, c(2, 3, 4)])

# PPCM
row.names(inputList_PPCM) <- NULL
# inputList_PPCM = inputList_PPCM[!is.na(inputList_PPCM[,1])&!is.na(inputList_PPCM[,2]),]
inputList_PPCM <- distinct(inputList_PPCM, ENSEMBL, .keep_all = TRUE)
inputList_PPCM <- distinct(inputList_PPCM, SYMBOL, .keep_all = TRUE)
inputList_PPCM[, c(3, 4)] <- (allgenes_PPCM[allgenes_PPCM$EnsemblID %in% inputList_PPCM[, 1], c(2, 3)])
inputList_PPCM <- data.frame(inputList_PPCM[, c(2, 3, 4)])

row.names(inputList_INTdcm) <- NULL
# inputList_INTdcm = inputList_INTdcm[!is.na(inputList_INTdcm[,1])&!is.na(inputList_INTdcm[,2]),]
inputList_INTdcm <- distinct(inputList_INTdcm, ENSEMBL, .keep_all = TRUE)
inputList_INTdcm <- distinct(inputList_INTdcm, SYMBOL, .keep_all = TRUE)
inputList_INTdcm[, c(3, 4)] <- (allgenes_INTdcm[allgenes_INTdcm$EnsemblID %in% inputList_INTdcm[, 1], c(2, 3)])
inputList_INTdcm <- data.frame(inputList_INTdcm[, c(2, 3, 4)])

row.names(inputList_INThcm) <- NULL
# inputList_INThcm = inputList_INThcm[!is.na(inputList_INThcm[,1])&!is.na(inputList_INThcm[,2]),]
inputList_INThcm <- distinct(inputList_INThcm, ENSEMBL, .keep_all = TRUE)
inputList_INThcm <- distinct(inputList_INThcm, SYMBOL, .keep_all = TRUE)
inputList_INThcm[, c(3, 4)] <- (allgenes_INThcm[allgenes_INThcm$EnsemblID %in% inputList_INThcm[, 1], c(2, 3)])
inputList_INThcm <- data.frame(inputList_INThcm[, c(2, 3, 4)])

#-------------------------------------------------------------------------------
# Creating expression matrix with gene symbols
#-------------------------------------------------------------------------------
exp_mat <- ABgxData
GeneSymbols <- bitr(row.names(exp_mat),
                    fromType = "ENSEMBL",
                    toType = c("SYMBOL"),
                    OrgDb = org.Hs.eg.db
)
row.names(GeneSymbols) <- NULL
GeneSymbols <- distinct(GeneSymbols, ENSEMBL, .keep_all = TRUE)
GeneSymbols <- distinct(GeneSymbols, SYMBOL, .keep_all = TRUE)
row.names(GeneSymbols) <- GeneSymbols[, 2]
GeneSymbols <- cbind(GeneSymbols, ABgxData[row.names(ABgxData) %in% GeneSymbols$ENSEMBL, ])
GeneSymbols <- GeneSymbols[, -c(1, 2)]

save(inputList_DCM,inputList_HCM,inputList_PPCM,inputList_INTdcm,
     inputList_INThcm,GeneSymbols, file = "ORA and enrichment analysis - Input Data.RData")
#-------------------------------------------------------------------------------
# 5. Pathway enrichment analysis 
#-------------------------------------------------------------------------------
# Loading saved data
load(file= "ORA and enrichment analysis - Input Data.RData")

#-------------------------------------------------------------------------------
# PathfindR 
#-------------------------------------------------------------------------------

dir <- ""
for (i in c("inputList_DCM", "inputList_HCM", "inputList_PPCM", "inputList_INTdcm", "inputList_INThcm")){
  if (i == "inputList_DCM") {
    dir <- "C:/Users/kaila/Documents/R/CardiomyopathyGroup2/pathfindR_results/DCM"
    outputVar <- "pathfindResult_DCM"
  } else {
    if (i == "inputList_HCM") {
      dir <- "C:/Users/kaila/Documents/R/CardiomyopathyGroup2/pathfindR_results/HCM"
      outputVar <- "pathfindResult_HCM"
    } else {
      if (i == "inputList_PPCM"){
      dir <- "C:/Users/kaila/Documents/R/CardiomyopathyGroup2/pathfindR_results/PPCM"
      outputVar <- "pathfindResult_PPCM"
      } else {
      if (i == "inputList_INTdcm"){
        dir <- "C:/Users/kaila/Documents/R/CardiomyopathyGroup2/pathfindR_results/INTdcm"
        outputVar <- "pathfindResult_intDCM"
      } else {
        dir <- "C:/Users/kaila/Documents/R/CardiomyopathyGroup2/pathfindR_results/INThcm"
        outputVar <- "pathfindResult_intHCM"
      }
      }
    }
  }
  assign(
    outputVar,
    run_pathfindR(get(i),
                  p_val_threshold = 0.05,
                  convert2alias = TRUE,
                  pin_name_path = "STRING", # Name of the chosen PIN - "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG"
                  sig_gene_thr = 0.02, # Threshold for minimum proportion of significant genes used in filtering active subnetworks
                  gene_sets = "KEGG", # gene set used for enrichment analysis - "KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC", "GO-MF"
                  min_gset_size = 10,
                  max_gset_size = 500,
                  output_dir = dir,
                  plot_enrichment_chart = TRUE
    )
  )
}

#-------------------------------------------------------------------------------
# PathfindR results 
#-------------------------------------------------------------------------------
#UpSet_plot(pathfindResult_DCM,num_terms = 10)
pdf("pathfindR_results.pdf")
term_gene_graph(pathfindResult_PPCM, num_terms = 5, node_size = "num_genes", use_description = T)
term_gene_graph(pathfindResult_DCM, num_terms = 5, node_size = "num_genes", use_description = T)
term_gene_graph(pathfindResult_HCM, num_terms = 5, node_size = "num_genes", use_description = T)
graphics.off()

# knitr::kable(head(pathfindResult_DCM, 2))

# downloads an illustration of the pathway
visualize_terms(
  result_df = pathfindResult_DCM,
  input_processed = input_processing(
    input = inputList_DCM, # the input: in this case, differential expression results
    p_val_threshold = 0.05, # p value threshold to filter significant genes
    pin_name_path = "KEGG", # the name of the PIN to use for active subnetwork search
    convert2alias = TRUE # boolean indicating whether or not to convert missing symbols to alias symbols in the PIN
  ),
  hsa_KEGG = TRUE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
  pin_name_path = "STRING"
)

# plotting enrichment charts
pdf("Enrichment charts.pdf")
enrichment_chart(pathfindResult_DCM, top_terms = 15)
enrichment_chart(pathfindResult_HCM, top_terms = 15)
enrichment_chart(pathfindResult_PPCM, top_terms = 15)
graphics.off()

# clustering enriched terms
diseases <- c("PPCM") #c("DCM", "HCM", "PPCM")
for (i in diseases) {
  outputVar <- paste(i, "_clustered", sep = "")
  source <- get(paste("pathfindResult_", i, sep = ""))
  assign(outputVar, cluster_enriched_terms(source, plot_hmap = F, plot_dend = F, plot_clusters_graph = F))
}

#knitr::kable(head(DCM_clustered, 2))
# For hsa_KEGG = TRUE,  KEGG  human  pathway  diagrams  are  created,  affected
# nodes  colored  by up/down regulation status. For other gene sets, interactions
# of affected genes are determined (via a shortest-path algorithm) and are
# visualized (colored by change status) using igraph.

# saves visualizations in the folder "term_visualizations" under the current working directory

## Vector of "Case" IDs
casesPPCM <- row.names(sampleInfo[sampleInfo$Disease == "PPCM", ])
casesDCM <- row.names(sampleInfo[sampleInfo$Disease == "DCM", ])
casesHCM <- row.names(sampleInfo[sampleInfo$Disease == "HCM", ])
healthy <- row.names(sampleInfo[sampleInfo$Disease == "Donor", ])

## Calculate scores for representative terms and plot heat map using term descriptions
for (i in c( "DCM", "HCM","PPCM")) {
  enrichment <- get(paste(i,"_clustered",sep=""))
  cases <- get(paste("cases",i,sep=""))
  outputVar <- paste("score_matrix_",i,sep="")
  assign(outputVar, score_terms(
    enrichment_table = enrichment[enrichment$Status == "Representative", ],
    exp_mat = as.matrix(GeneSymbols),
    cases = cases,
    control_title = "healthy", # default = "Control"
    use_description = TRUE, # default FALSE
    label_samples = FALSE, # default = TRUE
    case_title = as.character(i), # default = "Case"
    low = "#f7797d", # default = "green"
    mid = "#fffde4", # default = "black"
    high = "#1f4037") # default = "red"
  )
}

# saving the result data
save(pathfindResult_DCM,pathfindResult_HCM,pathfindResult_PPCM, DCM_clustered, HCM_clustered, PPCM_clustered, score_matrix_DCM, score_matrix_HCM, score_matrix_PPCM, file = "pathfindR results.RData")

#-------------------------------------------------------------------------------
# Creating data for WebGestalt - an online tool for enrichment analysis
#-------------------------------------------------------------------------------

# df of ordered FC with names (Ensembl) - for WebGestalt GSEA
dcmVectorGSEA <- allgenes_DCM[, c(1, 2)]
hcmVectorGSEA <- allgenes_HCM[, c(1, 2)]
ppcmVectorGSEA <- allgenes_PPCM[, c(1, 2)]

write.table(dcmVectorGSEA, file = "dcmVectorGSEA.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hcmVectorGSEA, file = "hcmVectorGSEA.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ppcmVectorGSEA, file = "ppcmVectorGSEA.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# list of ensembl names for DCM HCM PPCM INTdcm and INThcm DE genes - for ORA
dcmDEensembl <- DEGdcmEnsembl[, 1]
hcmDEensembl <- DEGhcmEnsembl[, 1]
ppcmDEensembl <- DEGppcmEnsembl[, 1]
INTdcmDEensembl <- DEGdcmINTEnsembl[, 1]
INThcmDEensembl <- DEGhcmINTEnsembl[, 1]

write(dcmDEensembl, file = "dcmtopDEensembl.txt")
write(hcmDEensembl, file = "hcmtopDEensembl.txt")
write(ppcmDEensembl, file = "ppcmtopDEensembl.txt")
write(row.names(gxData), file = "universeGenesEnsembl.txt")

save(dcmVectorGSEA,hcmVectorGSEA,ppcmVectorGSEA,dcmDEensembl,hcmDEensembl,ppcmDEensembl,INTdcmDEensembl,INThcmDEensembl,file="TopGO - input data.RData")
#-------------------------------------------------------------------------------
# topGO - for GO/KEGG enrichment analysis
#-------------------------------------------------------------------------------
# Load saved data
lnames <- load(file="TopGO - input data.RData")

# topGO annotation and enrichment analysis for all diseases across all GO ontologies. 
diseases <- c("DCM", "HCM", "PPCM")
for (i in diseases){
  InterestingGenesforTopGO <- get(paste(tolower(i),"DEensembl", sep="")) # Vector of interesting Ensembl IDs (such as DE genes)
  for (ont in c("BP", "CC", "MF")) {
    outputVar <- paste("allRes",i,ont,sep="_") #setting output variable name
    file <- paste("TopGO_Results_",i,".xlsx", sep="") #setting file name
    print(ont)
    dataForTopGO <- annFUN.org(ont, mapping = "org.Hs.eg.db", ID = "ensembl")
    universeGenes <- unique(unlist(dataForTopGO))
    factorGenesTopGO <- factor(as.integer(universeGenes %in% InterestingGenesforTopGO))
    names(factorGenesTopGO) <- universeGenes
    GOdata <- new("topGOdata",
                  ontology = ont, allGenes = factorGenesTopGO, nodeSize = 5,
                  annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl"
    )
    goFisherResult <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher") # recommended setting	= parentchild algorithm
    #goKSresults <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    allRes <- GenTable(GOdata,
                       Pvalue = goFisherResult, # classicFisher = goFisherResult,classicKS = goKSresults,
                       topNodes = length(goFisherResult@score),
                       orderBy = "Pvalue", ranksOf = "Pvalue",
                       numChar = 1E9
    )
    if (all(as.numeric(allRes$Pvalue) >= 0.05)) {
      allRes <- allRes[1, ]
      allRes[, 1:6] <- NA
    } else {
      allRes <- allRes[1:max(which(as.numeric(allRes[, 6]) < 0.05)), ] # leave only results with p < 0.05
    }
    
    write.xlsx(allRes,
               file = file,
               sheet = ifelse(ont == "BP", "BiologicalProcess",
                              ifelse(ont == "CC", "CellularComponent", "MolecularFunction")
               ),
               row.names = FALSE,
               append = TRUE
    )
  }
}


save(allRes_DCM_BP, allRes_DCM_CC, allRes_DCM_MF, allRes_HCM_BP, allRes_HCM_CC, allRes_HCM_MF, allRes_PPCM_BP, allRes_PPCM_CC, allRes_PPCM_MF, file = "TopGO result.RData")
#--------------------------#
#Visualisation 

#histogram of p values
hist(score(goFisherResult), 50, xlab = "p-values") 

# GO term structure plot. 
showSigOfNodes(GOdata, score(goFisherResult), firstSigNodes = 5, useInfo = "all")

# The subgraph induced by the top 5 GO terms identified by the parentchild algorithm for
# scoring GO terms for enrichment.  Rectangles indicate the 5 most significant terms.
# Rectangle color represents the relative significance,ranging from dark red
# (most significant) to bright yellow (least significant). For each node, some
# basic information is displayed.  The first two lines show the GO identifier and
# a trimmed GO name.  In the third line the raw p-value is shown.  The forth line
# is showing the number of significant genes and the total number of genes
# annotated to the respective GO term.

# Saving the GO term structure plot
printGraph(GOdata, goFisherResult, firstSigNodes = 10, fn.prefix = "Significant nodes ppcmDEensembl", useInfo = "all", pdfSW = T)


# Choosing GO term of interest to analyse individually
goID <- allRes[1, "GO.ID"] 

# Plots density curve to show distribution of gene scores of genes within a GO term. 
#print(showGroupDensity(GOdata, goID, ranks = TRUE))