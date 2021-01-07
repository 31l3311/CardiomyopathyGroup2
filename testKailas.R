# Importing required libraries and packages
require("ggplot2")
require("tidyverse")
require("tidyr")
require("pcaMethods")
require("Biobase")
require("limma")
require("biomaRt")
library(WGCNA)
options(stringsAsFactors = FALSE)

if (!requireNamespace("gridExtra", quietly = TRUE))
  install.packages("gridExtra")
require("gridExtra")

if (!requireNamespace("viridis", quietly = TRUE))
  install.packages("viridis")
require("viridis")

#-----------------------------------------------------------------------------#
# Assignment 1: Importing the data and inspecting sample information
#-----------------------------------------------------------------------------#

# Assignment 1a - Importing data

setwd("C:/Users/kaila/Documents/MSc/Period 2/Experimental methods and data management/R skill sessions/Data") # setting working directory
gxData <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", row.names = 1, header = TRUE) # importing gene expression data
sampleInfo <- read.delim("MAGNET_SampleData_19112020.txt", row.names = 1, header = TRUE) # importing phenotype data
geneTotExonLengths <- read.delim("MAGNET_exonLengths.txt", as.is = T, row.names = 1) # importing gene length data
gxT<- t(gxData) #Transpose of expression data required for WGCNA

gsg=goodSamplesGenes(gxT) #Checks any missing values

#---------------------------------------------------------------------------------
# Assignment 3a - differential gene expression analysis
#---------------------------------------------------------------------------------
# Creating an expressionSet object
eset <- ExpressionSet(as.matrix(gxData),
                      phenoData = AnnotatedDataFrame(sampleInfo),
                      featureData = AnnotatedDataFrame(geneTotExonLengths)
)

# Creating a design matrix with no intercept and 'Disease', 'Sex' and 'Ethnicity' 
# as explanatory variables. All covariates were included in the linear model to 
# account for any clustering that may be present in the dataset which may be 
# visible only in 5th or later principle components. However, only disease
# was seen to the reason for clustering in the first 5 principle components. 

design <- model.matrix(~ 0 + Disease+ Sex + Ethnicity,  data = pData(eset))

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
fit2 <- eBayes(fit2)

# Classifying genes as up, down or not significantly regulated.
results <- decideTests(fit2)
summary(results)
t4<-summary(results)

# Differentially expressed genes in DCM vs Donors- logFC threshold was kept to log2(1.5)
# and Benjamini-Hochberg correction was applied to p value. 0.05 was set as the cutoff
# for this adjusted p value (False Discovery Rate).

deg <- topTable(fit2,
                coef = "DCMvDon", adjust.method = "BH", lfc = log2(1.5),
                number = nrow(fit2), sort.by = "none", p.value = 0.05
)

# Subsetting up and down regulated genes.
upRegulatedGenes<-deg[deg$logFC>0,] 
downRegulatedGenes<-deg[deg$logFC<0,]

# Histogram of p values as a sanity check
hist(deg[, "adj.P.Val"], main="Adjusted P Value")

# volcano plot with top 5 DE genes labeled.
volcanoplot(fit2, highlight = 5, names = rownames(fit2$genes), main = "DCM vs Donor")


newdata=as.data.frame(t(gxData[which(results[,"DCMvDon"]!=0),]))
newTrait = sampleInfo[rownames(newdata),c(1,2,3)]
#---------------------------------------------------------------------------------------

#check if there are outliers by clustering the samples 
sampleTree = hclust(dist(gxT), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="",
     xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

gxTrait = sampleInfo

#all entries must be numeric for further analysis
gxTrait$Age <- as.numeric(gxTrait$Age)
gxTrait=gxTrait[,-4]
gxTrait[gxTrait=="Male"]=0
gxTrait[gxTrait=="Female"]=1
gxTrait[gxTrait=="Donor"]=0
gxTrait[gxTrait=="DCM"]=1
gxTrait[gxTrait=="HCM"]=2
gxTrait[gxTrait=="PPCM"]=3
gxTrait <- mutate_all(gxTrait, function(x) as.numeric(as.character(x)))
collectGarbage()

traitColours= numbers2colors(gxTrait,signed=F)

#white is a low value, red a high value. This shows relation bw sample info and dendrogram
#can see if there is any clustering of the sample info. Here there don't seem to be any significant
#ones other than that white patch in the right half of disease - which means 
#those are all donors. 

plotDendroAndColors(sampleTree, traitColours,
                    groupLabels = names(gxTrait), cex.dendroLabels = 0.5, 
                    main = "Sample dendrogram and trait heatmap")



powers = seq(1,15, by=1)
sft = pickSoftThreshold(gxT, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,
     signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net = blockwiseModules(gxT, maxBlockSize = 11000, power = 8,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "expTOM", 
                       verbose = 3)

# open a graphics window
sizeGrWindow(15, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.01,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs

# Define numbers of genes and samples
nGenes = ncol(gxT);
nSamples = nrow(gxT);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(gxT, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, gxTrait, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


sizeGrWindow(20,20)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep ="");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(gxTrait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Define variable time containing the time column of datTrait
disease = as.data.frame(gxTrait$Disease);
names(disease) = "disease"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(gxT, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(gxT, disease, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(disease), sep="");
names(GSPvalue) = paste("p.GS.", names(disease), sep="");


#yellow brown and grey are strongly associated as seen in the heatmap. so lets explore those

#module = "red"
module = "yellow"
#module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for disease",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# Create the starting data frame
geneInfo0 = data.frame(Gene.ID = colnames(newdata),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for time
modOrder = order(-abs(cor(MEs, disease, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.disease));
geneInfo = geneInfo0[geneOrder, ]

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(data.filtered.dcm, power = 6);
#save(TOM, file = "WGCNA-TOM.RData")

# Select modules
#modules = c("black");
#modules = c("brown");
modules = c("red");
# Select module probes
probes = names(data.filtered.dcm)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])


ensID = as.matrix(names(newdata))
colnames(ensID) <- c("ensembl_gene_id")
ensID2annotate = as.matrix(ensID[moduleColors=="yellow"])

intModules=c("brown", "yellow", "grey")
for (module in intModules){
  ensID2annotate = ensID[moduleColors==module]
  entrezID = getBM(attributes = c("ensembl_gene_id","entrezgene_id",
  ),
  filters = "ensembl_gene_id",
  values = ensID2annotate,
  mart = ensembl
  )
  fileName = paste("Entrez IDs-", module, ".txt", sep="");
  write.table(as.data.frame(entrezID), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

GOenr = GOenrichmentAnalysis(moduleColors, geneDetails[,3], organism = "human", nBestP = 10);


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(newdata, power = 8);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

ensID2annotate<-as.matrix(ensID)

# Setting ensembl mart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Getting the gene ID, name and description from ensembl database
geneDetails <- getBM(
  attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"), filters = "ensembl_gene_id",
  values = ensID2annotate,
  mart = ensembl
)
# Merging files. Gene IDs that do not have entries on Ensembl will have values = NA
geneDetails <- merge(geneIDs, geneDetails, by.x = "ensembl_gene_id", all.x = TRUE)

# Checking if row names match
all(rownames(geneDetails) == rownames(geneIDs))





