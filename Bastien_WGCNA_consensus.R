#######################################################################################
#======================================================================================
# 
# WGCNA analysis
# PRO4002:Project 1, group 2.
# Author: Bastien Nihant
# Date of last update: 14-01-2021
#
#======================================================================================
#######################################################################################


# Display the current working directory
getwd();
setwd("D:/Documents/Github/Research_project_1/CardiomyopathyGroup2")

# Load the package
library(WGCNA);
library(dplyr)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Enable multi-threading if the tool used for R supports it.
#enableWGCNAThreads()
#=====================================================================================
#
#  Opening data files
#
#=====================================================================================



#Read the entire data set
data <- read.delim("D:/Documents/Github/Research_project_1/ABgxData.txt", row.names = 1)
#read the entire metadata
traitData <- read.delim("D:/Documents/Github/Research_project_1/MAGNET_SampleData_19112020.txt", row.names = 1)
# Take a quick look at what is in the data sets (caution, longish output):
dim(data)
names(data)


#=====================================================================================
#
#  Setting up the expression data in a usable format
#
#=====================================================================================


# We work with three sets:
nSets = 3;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("DCM_Donor", "HCM_Donor", "PPCM_Donor")
shortLabels = c("DCM", "HCM", "PPCM")

#Divide the Donor group. 
#!!!!Due to randomization this part should not be run anymore if possible. Last ran 18th of January 2021 at 16:37.
nfemales <- sum(traitData$Sex == "Female"| traitData$Disease == "Donor")
donors <- traitData[traitData$Disease == "Donor",]
donorPPCM <- donors[donors$Sex == "Female",]
donorPPCM <- donorPPCM[sample(1:nrow(donorPPCM), 8),]
traitDataIntermediate <- donors[!(row.names(donors) %in% row.names(donorPPCM)),]
donorDCM <- traitDataIntermediate[sample(1:nrow(traitDataIntermediate), 115),]
donorHCM <- traitDataIntermediate[!(row.names(traitDataIntermediate) %in% row.names(donorDCM)),]
donorDCMIndex <- colnames(data) %in% row.names(donorDCM) 
donorHCMIndex <- colnames(data) %in% row.names(donorHCM) 
donorPPCMIndex <- colnames(data) %in% row.names(donorPPCM) 

save(donorDCMIndex,donorHCMIndex,donorPPCMIndex, file = "Intermediate_data/Donors_Index.RData")
# Form multi-set expression data: columns starting from 9 contain actual expression data
# multiExpr = vector(mode = "list", length = nSets)
# 
# multiExpr[[1]] = list(data = as.data.frame(t(data[,traitData$Disease == "Donor"|traitData$Disease == "DCM"])));
# rownames(multiExpr[[1]]$data) = colnames(data[,traitData$Disease == "Donor"|traitData$Disease == "DCM"]);
# multiExpr[[2]] = list(data = as.data.frame(t(data[,traitData$Disease == "Donor"|traitData$Disease == "HCM"])));
# rownames(multiExpr[[2]]$data) = colnames(data[,traitData$Disease == "Donor"|traitData$Disease == "HCM"]);
# multiExpr[[3]] = list(data = as.data.frame(t(data[,traitData$Disease == "Donor"|traitData$Disease == "PPCM"])));
# rownames(multiExpr[[3]]$data) = colnames(data[,traitData$Disease == "Donor"|traitData$Disease == "PPCM"]);
# 
# exprSize = checkSets(multiExpr)

indexes <- load(file = "Intermediate_data/Donors_Index.RData")
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(data[,donorDCMIndex|traitData$Disease == "DCM"])));
rownames(multiExpr[[1]]$data) = colnames(data[,donorDCMIndex|traitData$Disease == "DCM"]);
multiExpr[[2]] = list(data = as.data.frame(t(data[,donorHCMIndex|traitData$Disease == "HCM"])));
rownames(multiExpr[[2]]$data) = colnames(data[,donorHCMIndex|traitData$Disease == "HCM"]);
multiExpr[[3]] = list(data = as.data.frame(t(data[,donorPPCMIndex|traitData$Disease == "PPCM"])));
rownames(multiExpr[[3]]$data) = colnames(data[donorPPCMIndex|traitData$Disease == "PPCM"]);

exprSize = checkSets(multiExpr)


#=====================================================================================
#
#  Checking for for missing values, removing samples with excessive NAs
#
#=====================================================================================


# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK


if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}


#=====================================================================================
#
#  Creating the sample clustering tree for each set
#
#=====================================================================================


sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}



pdf(file = "Results/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();


#=====================================================================================
#
#  Removing outliers. This step was skipped.
#
#=====================================================================================

# 
# # Choose the "base" cut height for the first data set
# baseHeight = 16
# # Adjust the cut height for the male data set for the number of samples
# cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1]);
# # Re-plot the dendrograms including the cut lines
# pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
# par(mfrow=c(2,1))
# par(mar = c(0, 4, 2, 0))
# for (set in 1:nSets)
# {
#   plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
#        xlab="", sub="", cex = 0.7);
#   abline(h=cutHeights[set], col = "red");
# }
# dev.off();



# 
# for (set in 1:nSets)
# {
#   # Find clusters cut by the line
#   labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
#   # Keep the largest one (labeled by the number 1)
#   keep = (labels==1)
#   multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
# }
# collectGarbage();
# # Check the size of the leftover data
# exprSize = checkSets(multiExpr)
# exprSize


#=====================================================================================
#
#  Setting up the metadata in a usable format
#
#=====================================================================================


dim(traitData)
names(traitData)
datTraits<- traitData

datTraits <- traitData
datTraits$Age <- as.numeric(datTraits$Age)
datTraits <- datTraits[,-2]
datTraits[datTraits=="Donor"]<-0
datTraits[datTraits=="DCM"]<-1
datTraits[datTraits=="HCM"]<-1
datTraits[datTraits=="PPCM"]<-1
datTraits[datTraits=="Caucasian"]<-0
datTraits[datTraits=="African.American"]<-1
allTraits <- mutate_all(datTraits, function(x) as.numeric(as.character(x)))
row.names(allTraits) <- row.names(traitData)

collectGarbage();

# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, row.names(allTraits));
  Traits[[set]] = list(data = allTraits[traitRows, ]);
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;


#=====================================================================================
#
#  Saving the usable files to avoid having to re-run the previous part
#
#=====================================================================================


save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, 
     file = "Intermediate_data/Consensus-dataInput.RData");
#=====================================================================================
#
#  Loading data from the previous part
#
#=====================================================================================

# Load the data saved in the first part
lnames = load(file = "Intermediate_data/Consensus-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#=====================================================================================
#
# Choosing thresholding power.
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red", "blue")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "Results/scaleFreeAnalysisFinal.pdf", wi = 8, he = 6);
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();


#=====================================================================================
#
#  Running Blockwise consensus analysis.
#
#=====================================================================================


bnet = blockwiseConsensusModules(
  multiExpr, maxBlockSize = 8000, power = 9, minModuleSize = 30,
  deepSplit = 2, 
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveIndividualTOMs = TRUE, 
  saveConsensusTOMs = TRUE, verbose = 5)
save(bnet, file = "Intermediate_data/WGCNA-net.RData")
consMEs = bnet$multiMEs
moduleLabels = bnet$colors;
moduleColors = labels2colors(moduleLabels)
consTree = bnet$dendrograms[[1]];
#=====================================================================================
#
#  Loading data from the blockwise analysis.
#
#=====================================================================================
# 
# 
load(file = "Intermediate_data/WGCNA-net.RData")
bwLabels = matchLabels(bnet$colors, bnet$colors, pThreshold = 1e-7);
bwColors = labels2colors(bwLabels)

#=====================================================================================
#
# Examining the blockwise dendograms.
#
#=====================================================================================

sizeGrWindow(12,6)
pdf(file = "Results/BlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6);
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4);
nBlocks = length(bnet$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks)
  plotDendroAndColors(bnet$dendrograms[[block]], moduleColors[bnet$blockGenes[[block]]],
                      "Module colors",
                      main = paste("Gene dendrogram and module colors in block", block),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
dev.off();


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================
# 
# 
# sizeGrWindow(12,9);
# pdf(file="Plots/SingleDendro-BWColors.pdf", wi = 12, he = 9);
# plotDendroAndColors(consTree,
#                     cbind(moduleColors, bwColors),
#                     c("Single block", "Blockwise"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Single block consensus gene dendrogram and module colors")
# dev.off();
#=====================================================================================

#=====================================================================================
#
#  Calculating Module-trait correlations.
#
#=====================================================================================


# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}


#=====================================================================================
#
#  Plotting heatmaps of module-traits relationships
#
#=====================================================================================


# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,7)
pdf(file = "Results/ModuleTraitRelationships-DCM.pdf", wi = 10, he = 7);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                    signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();
# Plot the module-trait relationship table for set number 2
set = 2
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                    signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "Results/ModuleTraitRelationships-HCM.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();

# Plot the module-trait relationship table for set number 3
set = 3
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                    signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "Results/ModuleTraitRelationships-PPCM.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();

#=====================================================================================
#
#  Calculating consensus module-traits correlations.
#
#=====================================================================================


# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative], moduleTraitCor[[3]][negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative], moduleTraitPvalue[[3]][negative]);
# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive], moduleTraitCor[[3]][positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive], moduleTraitPvalue[[3]][positive]);


#=====================================================================================
#
#  Plotting consensus module-traits correlations.
#
#=====================================================================================


textMatrix =  paste(signif(consensusCor, 2), "\n(",
                    signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "Results/ModuleTraitRelationships-consensus1.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[1]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste0("Consensus module--trait relationships across"))
dev.off();                           



#=====================================================================================
#
#  Setting up result matrices
#
#=====================================================================================

# Retrieving the gene names
probes = row.names(data)



consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}



GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z + GS[[3]]$Z)/sqrt(2);
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z + kME[[3]]$Z)/sqrt(2);
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);



GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[3]]$cor, GS[[1]]$p, GS[[2]]$p, GS[[3]]$p, GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 8*nTraits)
rownames(GSmat) = rownames(data);
colnames(GSmat) = spaste(
  c("GS.set1.", "GS.set2." ,"GS.set3.", "p.GS.set1.", "p.GS.set2.", "p.GS.set3.", "Z.GS.meta.", "p.GS.meta"),
  rep(traitNames, rep(8, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[3]]$cor, kME[[1]]$p, kME[[2]]$p, kME[[3]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 8*nMEs)
rownames(kMEmat) = rownames(data);
colnames(kMEmat) = spaste(
  c("kME.set1.", "kME.set2.", "kME.set3.", "p.kME.set1.", "p.kME.set2.", "p.kME.set3.", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(8, nMEs)))


#=====================================================================================
#
#  Creating and saving final result table
#
#=====================================================================================


info = data.frame(ensemblID = rownames(data), ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  GSmat,
                  kMEmat);
write.table(info, file = "Results/consensusAnalysis-CombinedNetworkResults.txt", sep = "\t",
          row.names = FALSE, quote = FALSE)

##############
#
#Testing to see which module contains the most DEGs
#
#########
modulesOfInterests <- c("darkred", "red", "black","purple")
#Opening complete results table from DE analysis and keeping all DEGs
degRes <- read.table("results_all.txt", sep = "\t" )
degRes <- degRes[which(degRes$qvalue_DCM <= 0.05 | degRes$qvalue_HCM <= 0.05 | degRes$qvalue_PPC <= 0.05),]
#Testing for any DEG in midnightblue, magenta, royalblue, purple, salmon, and turquoise

nDEG <- percDEG <- c()
for (n in c(1:length(modulesOfInterests))){
  nDEG[n] <- sum(rownames(degRes) %in% rownames(info[which(info$ModuleColor == modulesOfInterests[n]),]))
  percDEG[n] <- round(nDEG[n] / nrow(info[which(info$ModuleColor == modulesOfInterests[n]),]) * 100, digits = 1)
  print(paste0("Number of DEGs in ", modulesOfInterests[n] ," module (percentage) ",nDEG[n],"(",percDEG[n],")"))
}
#####
#Doing the same thing here but with genes that are DEGs for all diseases
####

#Opening complete results table from DE analysis and keeping overlap DEGs
degRes <- read.table("results_all.txt", sep = "\t" )
degRes <- degRes[which(degRes$qvalue_DCM <= 0.05 & degRes$qvalue_HCM <= 0.05 & degRes$qvalue_PPC <= 0.05),]
#Testing for any DEG in midnightblue, magenta, royalblue, purple, salmon, and turquoise

nDEG <- percDEG <- c()
for (n in c(1:length(modulesOfInterests))){
  nDEG[n] <- sum(rownames(degRes) %in% rownames(info[which(info$ModuleColor == modulesOfInterests[n]),]))
  percDEG[n] <- round(nDEG[n] / nrow(info[which(info$ModuleColor == modulesOfInterests[n]),]) * 100, digits = 1)
  print(paste0("Number of DEGs in ", modulesOfInterests[n] ," module (percentage) ",nDEG[n],"(",percDEG[n],")"))
}

save(list = ls(), file = "Results/environment_all.Rdata")
