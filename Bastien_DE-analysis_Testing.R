#------------------------------------#
#-----Loading required libraries-----#
#------------------------------------#

require("ggplot2")
require("tidyverse")
require("pcaMethods")
require("limma")
require("biomaRt")

#------------------------#
#-----Importing Data-----#
#------------------------#

# Fill in the directory in which any output should be saved.
outDir <- file.path("D:", "Documents", "Github", "Research_project_1", "Output")

# Importing gene expression data, sample information, and exon length. Make sure the "File_names" correspond to those in the code.
gxData <- read.delim("D:/Documents/Github/Research_project_1/MAGNET_GeneExpressionData_CPM_19112020.txt", row.names = 1)
sampleInfo <- read.delim("D:/Documents/Github/Research_project_1/MAGNET_SampleData_19112020.txt", row.names = 1, stringsAsFactors = TRUE)
geneTotExonLengths <- read.delim("MAGNET_exonLengths.txt", as.is = T, row.names = 1)

# Checking that the rows of sample_Info & geneTotExonLengths correspond to the columns of gxData.
all(row.names(sampleInfo) == colnames(gxData))
all(rownames(geneTotExonLengths) == rownames(gxData))

# Summarizing sampleInfo.
str(sampleInfo)
table(sampleInfo$Disease, sampleInfo$Sex)
table(sampleInfo$Disease, sampleInfo$Ethnicity)
ggplot(sampleInfo, aes(Age)) +
  geom_histogram(bins = 25) +
  facet_wrap(vars(Disease)) +
  theme_minimal()

# Transforming data into FPKM values.
cpm2fpkm <- function(x) {
  .t <- 2^(x) * 1E3 / geneTotExonLengths[, 1]
}
gxData_fpkm <- cpm2fpkm(gxData)


#--------------------------#
#-----Diagnostic plots-----#
#--------------------------#

# Making and saving a boxplot of gene expression values for each sample in the chosen disease groups (Donor and DCM by default).
# If the group contains more than 50 samples, it will be divided to produce readable graphs.
groupsOfInterest <- c("Donor", "DCM", "HCM", "PPCM")

for (group in groupsOfInterest) {
  plotData <- gather(gxData[, sampleInfo$Disease == group], key = "SampleID", value = "CPM")
  lGroup <- ncol(gxData[, sampleInfo$Disease == group])
  if (lGroup <= 50) {
    groupBoxplot <- ggplot(plotData, aes(SampleID, CPM)) +
      geom_boxplot() +
      coord_flip() +
      theme_classic() +
      ggtitle(group)
    print(groupBoxplot)
    ggsave(filename = paste0("gx_boxplot", "_", group, ".jpg"), plot = groupBoxplot, path = outDir)
  } else {
    nGroups <- (lGroup %/% 50) + 1 # This makes sure that no groups will be larger than 50
    groupSize <- ceiling(lGroup / nGroups)
    groups <- split(c(1:lGroup), ceiling(seq_along(c(1:lGroup)) / groupSize))
    for (g in c(1:nGroups)) {
      plotData2 <- gather(gxData[, sampleInfo$Disease == group][, unlist(groups[g])], key = "SampleID", value = "CPM")
      groupBoxplot <- ggplot(plotData2, aes(SampleID, CPM)) +
        geom_boxplot() +
        coord_flip() +
        theme_classic() +
        ggtitle(group)
      print(groupBoxplot)
      ggsave(filename = paste0("gx_boxplot", "_", group, g, ".jpg"), plot = groupBoxplot, path = outDir)
    }
  }
}

# Making and saving pca plot colored by each available covariates
pcaRes <- pca(t(gxData), nPcs = 10)
plotData <- cbind(data.frame(pcaRes@scores), sampleInfo)
for (g in colnames(sampleInfo)) {
  if (is.numeric(sampleInfo[, g]) == TRUE) {
    pcaGroupPlot <- ggplot(plotData, aes_string("PC7", "PC8", color = g)) +
      geom_point() +
      scale_color_gradient(low = "blue", high = "red")
  } else {
    pcaGroupPlot <- ggplot(plotData, aes_string("PC7", "PC8", color = g)) +
      geom_point()
  }
  print(pcaGroupPlot)
  ggsave(filename = paste0("PCA_plot3_", g, ".png"), plot = pcaGroupPlot, path = outDir)
}

#------------------------------#
#-----Statistical Analysis-----#
#------------------------------#

# Choosing the groups to be compared in the differential expression analysis (Donor and DCM by default) and subsetting sampleInfo and gxData accordingly.
group1 <- "Donor"
group2 <- "DCM"

#subSampleInfo <- sampleInfo[which(sampleInfo$Disease == group1 | sampleInfo$Disease == group2), ]
#subGxData <- gxData[, which(sampleInfo$Disease == group1 | sampleInfo$Disease == group2)]


# creating the design and contrast matrices.


# Note: based on the constructed PCA plots, ethnicity was determined to be a possible confounder to be included in the design.
# Other variables of sampelInfo did not have a visible impact according to the PCA plots.
# A different confounding variable can be selected by subsetting subSampleData with a different variable in the next line.
 ethnicity_ <- droplevels(sampleInfo$Ethnicity)
# names(confounder1) <- row.names(sampleInfo)
# sex_ <- droplevels(sampleInfo$Sex)
# names(confounder2) <- row.names(sampleInfo)
disease.sex <- paste(sampleInfo$Disease, sampleInfo$Sex, sep = ".")
design <- model.matrix(~ 0 + disease.sex + ethnicity_  )
colSums(design)

cm <- makeContrasts(
  DCM_Donor = (disease.sexDCM.Female + disease.sexDCM.Male) - (disease.sexDonor.Female + disease.sexDonor.Male),
  HCM_Donor = (disease.sexHCM.Female + disease.sexHCM.Male) - (disease.sexDonor.Female + disease.sexDonor.Male), 
  PPCM_Donor = disease.sexPPCM.Female - disease.sexDonor.Female, 
  InteractionDCM = (disease.sexDCM.Female - disease.sexDCM.Male) - (disease.sexDonor.Female - disease.sexDonor.Male),
  InteractionHCM = (disease.sexHCM.Female - disease.sexHCM.Male) - (disease.sexDonor.Female - disease.sexDonor.Male),
  levels = design
)
fit <- lmFit(gxData, design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2, trend = TRUE)

resultsSimple <- decideTests(fit2)
resultsDCM <- topTable(fit2, coef = 1, number = nrow(gxData), adjust.method = "BH", p.value = 0.05)
resultsHCM <- topTable(fit2, coef = 2, number = nrow(gxData), adjust.method = "BH", p.value = 0.05)
resultsPPCM <- topTable(fit2, coef = 3, number = nrow(gxData), adjust.method = "BH", p.value = 0.05)
resultsInteractionDCM <- topTable(fit2, coef = 4, number = nrow(gxData), adjust.method = "BH", p.value = 0.05)
resultsInteractionHCM <- topTable(fit2, coef = 5, number = nrow(gxData), adjust.method = "BH", p.value = 0.05)
# Outputting the number of up and downregulated genes, as well as the number of non-significant genes.
paste("DCM/  Up: ", nrow(resultsDCM[which(resultsDCM$logFC > 0), ]), ";Down: ", nrow(resultsDCM[which(resultsDCM$logFC < 0), ]), ";non-sig: ", nrow(gxData) - nrow(resultsDCM))
paste("HCM/  Up: ", nrow(resultsHCM[which(resultsHCM$logFC > 0), ]), ";Down: ", nrow(resultsHCM[which(resultsHCM$logFC < 0), ]), ";non-sig: ", nrow(gxData) - nrow(resultsHCM))
paste("PPCM/  Up: ", nrow(resultsPPCM[which(resultsPPCM$logFC > 0), ]), ";Down: ", nrow(resultsPPCM[which(resultsPPCM$logFC < 0), ]), ";non-sig: ", nrow(gxData) - nrow(resultsPPCM))
paste("InteractionDCM/  Up: ", nrow(resultsInteractionDCM[which(resultsInteractionDCM$logFC > 0), ]), ";Down: ", nrow(resultsInteractionDCM[which(resultsInteractionDCM$logFC < 0), ]), ";non-sig: ", nrow(gxData) - nrow(resultsInteractionDCM))
paste("InteractionHCM/  Up: ", nrow(resultsInteractionHCM[which(resultsInteractionHCM$logFC > 0), ]), ";Down: ", nrow(resultsInteractionHCM[which(resultsInteractionHCM$logFC < 0), ]), ";non-sig: ", nrow(gxData) - nrow(resultsInteractionHCM))
vennDiagram(resultsSimple[,c(1:3)])

degList <- row.names(gxData) %in% row.names(resultsDCM) | row.names(gxData) %in% row.names(resultsHCM) | row.names(gxData) %in% row.names(resultsPPCM)
gxDataDEG <- gxData[degList,]
write.table(gxDataDEG, file = "test_data.txt", sep = "\t", quote = FALSE, row.names = TRUE)

degOverlap <- rownames(gxData[row.names(gxData) %in% row.names(resultsDCM) & row.names(gxData) %in% row.names(resultsHCM) & row.names(gxData) %in% row.names(resultsPPCM),])
write(degOverlap, file = "DEG_overlap.txt")
degOverlapDCMHCM <- rownames(gxData[row.names(gxData) %in% row.names(resultsDCM) & row.names(gxData) %in% row.names(resultsHCM),])
write(degOverlapDCMHCM, file = "DEG_overlapDCMHCM.txt")

#making a table of results for all comparaisons
#First We need to get all the p.values and adjust for multiple comparaisons
pVal <- fit2$p.value
qval <- apply(pVal, 2, p.adjust, method = "BH")

#A dummy topTable is then called to get logFC values even for non-signifficant genes
logFC <- list()
for (n in c(1:5)){
  dummyTable <- topTable(fit2, coef = n, number = nrow(gxData), adjust.method = "BH", p.value = 1)
  logFC[[n]] <- dummyTable$logFC
}
remove(dummyTable)

resultsAll <- data.frame(logFC[[1]],pVal[,1],qval[,1],logFC[[2]],pVal[,2],qval[,2],
                         logFC[[3]],pVal[,3],qval[,3],logFC[[4]],pVal[,4],qval[,4],
                         logFC[[5]],pVal[,5],qval[,5])
colnames(resultsAll) <- c("logFC_DCM","pvalue_DCM","qvalue_DCM","logFC_HCM","pvalue_HCM","qvalue_HCM",
                          "logFC_PPCM","pvalue_PPCM","qvalue_PPC","logFC_SexDCM","pvalue_SexDCM","qvalue_SexDCM",
                          "logFC_SexHCM","pvalue_SexHCM","qvalue_SexHCM")
write.table(resultsAll, file = "results_all.txt", row.names = TRUE, quote = FALSE, sep ="\t")
#--------------------#
#-----Annotation-----#
#--------------------#

# Note: Genes which were not differentially expressed are not included in any further analysis, as they are considered "irrelevant".

# Cleaning the result table
results2 <- results[, c(1, 5)]

# Biomart query:
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name", "description"),
  values = row.names(results2),
  mart = ensembl
)
results2 <- merge.data.frame(results2, annotations, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)

# Cleaning the result table
row.names(results2) <- results2$Row.names
results2 <- results2[, -1]

#----------------------------------------------------#
#-----Relative gene expression / Noise estimation----#
#----------------------------------------------------#

# Note on the method: In order to identify noise, the FPKM values are taken to find in every female sample the median Y-chromosome gene expression.
# The gene whose expression is closest to this median in the largest number of samples is then selected as median noise gene (MNG).
# The mean expression of the MNG across samples is then calculated (in cpm) and called avgCpmNoise.
# Finally, the mean expression of every DEG is calculated (in cpm) and compared with avgCpmNoise.
# DEGs with an average expression that is smaller or larger than avgCpmNoise are labelled accordingly.

# Getting Y-chromosome genes
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

avgCpmNoise <- apply(gxData[medianNoiseGene, which(sampleInfo$Sex == "Female")], 1, FUN = mean)

# Calculating the mean expression of all genes and testing for noise
meanCpm <- apply(gxData, 1, FUN = mean)
noiseTestAll <- meanCpm > avgCpmNoise

#Removing too low genes based on noise
gxDataClean <- gxData[noiseTestAll,]

#flagging genes based on CEFIC rules -> no genes flagged so nice data
matCounts <- as.matrix(gxData)
diseases <- levels(sampleInfo$Disease)
threshold <- 1
qualityTest <- NULL
for (gene in row.names(gxData)){
  Check3<-c()
  for (s in diseases){
    Samples <- colnames(gxData[,which(sampleInfo[,3] == s)])
    matCounts2 <- matCounts[gene,Samples]
    Check1 <- sum(matCounts2 >= threshold) >= 0.75*length(Samples) 
    one <- max(matCounts2) - median(matCounts2)
    two <- sum(matCounts2)/(length(Samples)-1)
    Check2 <- one < two
    Check3 <- c(Check3, Check1 & Check2)
  }
  qualityTest <- c(qualityTest, sum(Check3)>=1)
}

#barplot
gene_to_plot <- "ENSG00000000003"
gene_count <- as.matrix(subset(gxData, rownames(gxData) == gene_to_plot))
barplot(gene_count, las=2, cex.names = 0.6)

# Calculating the mean expression of DEGs and testing for noise
meanCpm <- apply(gxData[row.names(resultsSimple), ], 1, FUN = mean)
noiseTest <- meanCpm > avgCpmNoise

#-------------------------------#
#-----Exporting the results-----#
#-------------------------------#

# Getting the mean expression of DEGs for both groups in cpm and FPKM
meanCmpDonor <- apply(gxData[row.names(results2), which(sampleInfo$Disease == group1)], 1, FUN = mean)
meanCmpDisease <- apply(gxData[row.names(results2), which(sampleInfo$Disease == group2)], 1, FUN = mean)

# Assembling the results in one table
resultsFinal <- data.frame(
  rownames(results2), results2$hgnc_symbol, results2$external_gene_name, results2$description, meanCmpDonor,
  meanCmpDisease, results2$logFC, results2$adj.P.Val, noiseTest
)

colnames(resultsFinal) <- c(
  "ensembl_gene_id", "HGNC_symbol", "external_gene_name", "description",
  "mean_CPM_Donor", "mean_CPM_Disease", "log_Fold_change", "adjusted_P-value",
  "higher_than_noise", "mean_FPKM_Donor", "mean_FPKM_Disease"
)

# Writing result table
write.table(resultsFinal, file = paste0(outDir, "/DE_Analysis_Results.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# Writing list of all measured genes (e.g. to use as background for over-representation analysis)
write(row.names(gxData), file = paste0(outDir, "/background_genes.txt"))
