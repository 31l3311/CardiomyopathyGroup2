
# Set up environment

#clear workspace and set string as factors to false
rm(list=ls())
options(stringsAsFactors = F)



# Install required packages


library(WGCNA)
library(rstudioapi)
library(dplyr)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# read data and combine them to input file for WGCNA

data <- read.table(file="test_data.txt", sep ="\t", header=TRUE)


# check if there are samples with missing data

gsg = goodSamplesGenes(data, verbose = 3);
gsg$allOK

# normalized counts from RNA-seq data should be log-transformed
data.log <- log2(data+1)
data.log <- as.data.frame(t(as.matrix(data.log)))
#cheching for outliers
sampleTree = hclust(dist(data.log), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
#Based on this we chose not to remove any outlier

traitData = read.table("MAGNET_SampleData_19112020.txt", header=TRUE, row.names = 1);

data.filtered.dcm = data.log[row.names(data.log) %in% row.names(traitData),]

# Form a data frame analogous to expression data that will hold the clinical traits.
# samples = rownames(data.filtered.dcm);
# traitRows = match(samples, traitData$SampleID);
# datTraits = traitData[traitRows, -1];
# rownames(datTraits) = traitData[traitRows, 1];

datTraits <- traitData
datTraits$Age <- as.numeric(datTraits$Age)
datTraits[datTraits=="Male"]<-0
datTraits[datTraits=="Female"]<-1
datTraits[datTraits=="Donor"]<-0
datTraits[datTraits=="DCM"]<-1
datTraits[datTraits=="HCM"]<-2
datTraits[datTraits=="PPCM"]<-3
datTraits[datTraits=="Caucasian"]<-0
datTraits[datTraits=="African.American"]<-1
datTraits <- mutate_all(datTraits, function(x) as.numeric(as.character(x)))

collectGarbage();


#merge together the filtered table with the information from the Giannakis dataset

# Cluster samples
sampleTree = hclust(dist(data.filtered.dcm), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
sizeGrWindow(12,12)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits), cex.dendroLabels = 0.5, 
                    main = "Sample dendrogram and trait heatmap")



save(data.filtered.dcm, datTraits, file = "WGCNA-input.RData")






#########################################
#Network construction and module detection
#########################################


# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "WGCNA-input.RData");
#The variable lnames contains the names of loaded variables.
lnames



# Choose a set of soft-thresholding powers
powers = seq(1,15, by=2)

# Call the network topology analysis function
sft = pickSoftThreshold(data.filtered.dcm, powerVector = powers, verbose = 5)

save(sft, file = "WGCNA-sft.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


tic <- Sys.time()
# looking at both - soft threshold and mean connectivity 
# I decided to go with power 6 for this small example dataset
net = blockwiseModules(data.filtered.dcm, power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "expTOM", 
                       verbose = 3)
tac <- Sys.time()
tac-tic
save(net, file = "WGCNA-net.RData")





# open a graphics window
sizeGrWindow(15, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "network-reconstruction.RData")




##########################################
Relate modules to external clinical traits
##########################################


# Load the expression and trait data saved in the first part
lnames = load(file = "WGCNA-input.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "network-reconstruction.RData");
lnames



# Define numbers of genes and samples
nGenes = ncol(data.filtered.dcm);
nSamples = nrow(data.filtered.dcm);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(data.filtered.dcm, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);



sizeGrWindow(20,20)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep ="");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
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
disease = as.data.frame(datTraits$Disease);
names(disease) = "disease"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(data.filtered.dcm, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(data.filtered.dcm, disease, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(disease), sep="");
names(GSPvalue) = paste("p.GS.", names(disease), sep="");



#module = "red"
module = "turquoise"
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
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,
                   abline = TRUE)



# Create the starting data frame
geneInfo0 = data.frame(Gene.ID = colnames(data.filtered.dcm),
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



write.csv(geneInfo, file = "geneInfo.csv", row.names = FALSE)



###################################
Network visualization 
###################################


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



