require("dplyr")
require("limma")
require("biomaRt")
require("xlsx")
require("ggplot2")
require("enrichplot")
require("clusterProfiler")
require("pathfindR")
require("pathview")
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
Ensdb <- EnsDb.Hsapiens.v86
options(stringsAsFactors = FALSE)


setwd("C:/Users/kaila/Documents/MSc/Period 2/Experimental methods and data management/R skill sessions/Data") 
#load(file = "usethis2401.RData")
#-------------------------------------------------------------------------------#
# 1 - Importing the data and inspecting sample information
#-------------------------------------------------------------------------------
# Importing data
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

# making file with all geneIDs
#write(row.names(gxData),file="backgroundGenes.txt")
# setting ensembl database
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#-------------------------------------------------------------------------------#
# 2 - Checking background noise level.
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

#-------------------------------------------------------------------------------#
# 3 - Annotating genes with entrez IDs, symbols and biotype
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

#creating gxData with rownames as symbols
symbols <- bitr(row.names(ABgxData),
                fromType = "ENSEMBL",
                toType = "SYMBOL",
                OrgDb = org.Hs.eg.db)
symbols <- distinct(symbols, ENSEMBL, .keep_all=TRUE)
symbols <- distinct(symbols, SYMBOL, .keep_all=TRUE)
row.names(symbols) <- symbols[,1]
ABgxSymbols <- ABgxData[row.names(symbols),]
row.names(ABgxSymbols) <- symbols[,2]

#biotype <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
#                 filters = "ensembl_gene_id", 
#                 values = row.names(ABgxData),
#                 mart = ensembl)
#write.table(ABgxData, file = "ABgxData.txt", sep="\t")
#save(sampleInfo, YGenes, ABgxData, ABgeneTotExonLengths, EntrezIDs, file = "DE analysis - Input Data.RData")
#-------------------------------------------------------------------------------#
# 4 - Differential gene expression analysis
#-------------------------------------------------------------------------------
# Load saved data 
lnames <- load(file = "DE analysis - Input Data.RData")

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
#deg_overlap <- row.names(ABgxData[row.names(ABgxData) %in% row.names(resDCM) & row.names(ABgxData) %in% row.names(resHCM) & row.names(ABgxData) %in% row.names(resPPCM), ])

# subsetting DE genes by disease
#overlapDEgenes <- ABgxData[deg_overlap, ]
# ppcmDEgenes = ABgxData[row.names(resPPCM),]
# dcmDEgenes = ABgxData[row.names(resDCM),]
# hcmDEgenes = ABgxData[row.names(resHCM),]

# Listing all significantly differentially expressed genes. 
# degList <- row.names(ABgxData) %in% row.names(resDCM) | row.names(ABgxData) %in% row.names(resHCM) | row.names(ABgxData) %in% row.names(resPPCM)
# ABgxDataDEG <- ABgxData[degList, ]

# Creating .txt files of allgenes_ dataframes for any online GSEA tools that may be used
#write.table(row.names(resDCM),file = "DEgenes_DCM.txt", sep="\t",quote=F, row.names = F)
#write.table(row.names(resHCM),file = "DEgenes_HCM.txt", sep="\t",quote=F, row.names = F)
#write.table(row.names(resPPCM),file = "DEgenes_PPCM.txt", sep="\t",quote=F, row.names = F)
#write.table(row.names(resIntDCM),file = "DEgenes_INTdcm.txt", sep="\t",quote=F, row.names = F)
#write.table(row.names(resIntHCM),file = "DEgenes_INThcm.txt", sep="\t",quote=F, row.names = F)
#-------------------------------------------------------------------------------#
# 5 - Logical vector to identify differentially expressed genes
#-------------------------------------------------------------------------------
listToDF <- function(list){
  sapply(list, "length<-", max(lengths(list)))
}

DEgenes <- list(row.names(resDCM), row.names(resHCM), row.names(resPPCM), row.names(resIntDCM), row.names(resIntHCM))
names(DEgenes) <- c("DEgenes_DCM", "DEgenes_HCM", "DEgenes_PPCM", "DEgenes_IntDCM","DEgenes_IntHCM")

# Create a data frame with DEGs, one column per disease
DEgenes<-listToDF(DEgenes)

DCM<-(row.names(ABgxData) %in% row.names(resDCM))
names(DCM) <- ABgxData
HCM<-(row.names(ABgxData) %in% row.names(resHCM))
names(HCM) <- ABgxData
PPCM<-(row.names(ABgxData) %in% row.names(resPPCM))
names(PPCM) <- ABgxData
DEgenesLogical <- data.frame(DCM,HCM,PPCM)
row.names(DEgenesLogical)<- row.names(ABgxData)
#DEgenesLogical["ENSG00000182866",]
# checking if the list made above contains all DE genes
#all((row.names(resDCM) %in% row.names(DEgenesLogical)[DEgenesLogical$DCM==TRUE]) &
#     row.names(DEgenesLogical)[DEgenesLogical$DCM==TRUE] %in% row.names(resDCM)) 

#write.table(DEgenesLogical, file="DEgenesLogical.txt", sep="\t", quote=FALSE)
#write.table(DEgenes, file="allDEgenes.txt", sep='\t', row.names = FALSE)

#-------------------------------------------------------------------------------#
# 6 - Creating dataframe for each disease of all gene IDs(ensembl), logfc, adj.p.value
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
#write.table(allgenes_DCM,file = "allgenes_DCM.txt", sep="\t")
#write.table(allgenes_HCM,file = "allgenes_HCM.txt", sep="\t")
#write.table(allgenes_PPCM,file = "allgenes_PPCM.txt", sep="\t")
#write.table(allgenes_INTdcm,file = "allgenes_INTdcm.txt", sep="\t")
#write.table(allgenes_INThcm,file = "allgenes_INThcm.txt", sep="\t")
#-------------------------------------------------------------------------------#
# 7 - GO and KEGG enrichment of WGCNA modules - clusterprofiler
#-------------------------------------------------------------------------------
# loading WGCNA results
WGCNAresults <- read.delim("consensusAnalysis-CombinedNetworkResults.txt", row.names = 1, header=TRUE)

# listing module colors 
modules <- unique(WGCNAresults[,2])

# listing modules of interest
modulesOfInterest <- c("tan", "darkred")#, "lightyellow","red","grey60","royalblue")

# Diseases of interest
diseases <- c("DCM", "HCM", "PPCM")

# Defining function to convert ENSEMBL to ENTREZID for each module
convertToEntrez <- function (inputlist){
  return(distinct(bitr(inputlist,"ENSEMBL","ENTREZID",OrgDb = org.Hs.eg.db),ENSEMBL,.keep_all=T))
}

# Defining function to convert ENSEMBL to SYMBOL for each module
convertToSymbol <- function(ensemblGeneList){
  return(distinct(bitr(ensemblGeneList,"ENSEMBL","SYMBOL",OrgDb = org.Hs.eg.db),ENSEMBL,.keep_all=T))
}

# listing genes in each module - can be used for ORA
for (module in modulesOfInterest){
  outputVar = paste0("genesIn_",module)
  assign(outputVar, row.names(WGCNAresults[WGCNAresults$ModuleColor==module,]))
  write(get(outputVar),file = paste0("genesIn_",module,".txt"))
}

# making dataframe of geneID, logFC, adjpvalue for each module for each disease 

for (module in modulesOfInterest){
  genesInModule <- get(paste0("genesIn_",module))
  for (disease in diseases){
    outputVar <- paste("fcData",module,disease,sep="_")
    source <- get(paste0("allgenes_",disease))
    assign(outputVar,data.frame(source[source$EnsemblID %in% genesInModule,]))
    dummy1= get(outputVar)
    row.names(dummy1)=dummy1[,1]
    assign(outputVar,dummy1)
  }
}
#write.table(fcData_darkred_DCM[,c(1,2)],file = "darkredDCMfc.txt",row.names=F,quote=F, col.names = F,sep=",")


# making named and ranked gene list of each module for each disease - for gsea on constituents of each module
for (module in modulesOfInterest){
  for (disease in diseases){
    source = get(paste("fcData",module,disease,sep="_"))
    outputVar = paste("rankedGeneList",module,disease,sep="_")
    assign(outputVar, source[,2])
    dummy1 = get(outputVar)
    names(dummy1) = row.names(source)
    dummy1 = sort(dummy1,decreasing = TRUE)
    assign(outputVar,dummy1)
  }
}

# making ENTREZID named and ranked gene list of each module for each disease - for gsea with KEGG
for (module in modulesOfInterest){
  for (disease in diseases){
    entNames <- convertToEntrez(get(paste0("genesIn_",module)))
    source <- get(paste("fcData",module,disease,sep="_"))
    column <- paste0(disease,"fc")
    dummy1<- source[source[,1] %in% entNames[,1],2]
    names(dummy1) <- entNames[,2]
    source <- sort((dummy1),decreasing=TRUE)
    outputVar <- paste("entrezRanked",module,disease,sep="_")
    assign(outputVar,source)
  }
}

# GO enrichment 
ontology <- c("BP", "CC", "MF", "ALL")
for (module in modulesOfInterest){
  source = get(paste0("genesIn_",module))
  file = paste0("goEnrichmentResults_",module,".xlsx")
  for (ont in ontology){
    dataForTopGO <- annFUN.org(ont, mapping = "org.Hs.eg.db", ID = "ensembl")
    universeGenes <- unique(unlist(dataForTopGO))
    outputVar = paste("ego",module,ont,sep="_")
    assign(outputVar, enrichGO(gene = source, 
                               universe = row.names(gxData),
                               keyType = "ENSEMBL",
                               OrgDb = org.Hs.eg.db,
                               ont = ont,
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.5,
                               qvalueCutoff = 0.5,
                               readable = FALSE,
                               minGSSize = 2,
                               maxGSSize = 500))
    resGO <- get(outputVar)
    write.xlsx(resGO,
               file = file,
               sheet = ifelse(ont == "BP", "BiologicalProcess",
                              ifelse(ont == "CC", "CellularComponent",
                                     ifelse(ont=="MF", "MolecularFunction","All"))
               ),
               row.names = FALSE,
               append = TRUE)
  }
}

# Kegg enrichment 
for (module in modulesOfInterest){
  source = get(paste0("genesIn_",module))
  dummy1 = bitr(source,"ENSEMBL", "ENTREZID", OrgDb = org.Hs.eg.db)
  dummy2 = distinct(dummy1,ENSEMBL,.keep_all = TRUE)
  source = dummy2[,2]
  file = paste0("KEGG_EnrichmentResults_",module,".xlsx")
  outputVar = paste("eKEGG",module,sep="_")
  assign(outputVar, enrichKEGG(gene = source, 
                               organism = "hsa", 
                               keyType = "ncbi-geneid",
                               pvalueCutoff = 0.5,
                               pAdjustMethod = "BH",
                               universe = EntrezIDs[,2],
                               minGSSize = 2,
                               maxGSSize = 500,
                               qvalueCutoff = 0.5))
    resKEGG <- get(outputVar)
    write.xlsx(resKEGG,
               file = file,
               row.names = FALSE)
}


#-------------------------------------------------------------------------------#
# 8 - Visualizing enrichment results
#-------------------------------------------------------------------------------
#barplots of GO enrichment 

pdf(file = paste0("Barplots of GO enrichment-",module,".pdf"))
for (module in modulesOfInterest){
  for (ont in ontologies){
    plotData = get(paste("ego",module,ont,sep="_"))
    a=barplot(plotData, drop=TRUE, showCategory = 10, 
            title = ifelse(ont=="BP", paste0("Biological Pathway-",module),
                          ifelse(ont=="CC",paste0("Cellular Component-",module),
                                 paste0("Molecular Function-",module))),font.size=8)+
      ggplot2::theme(legend.position =  "none",
                     plot.title = element_text(face="bold",hjust=0.5,vjust=0.5))
    print(a)
  }
}
dev.off()

#barplots of kegg enrichment
pdf(file = paste0("Barplots of KEGG enrichment-",module,".pdf"))
for (module in modulesOfInterest){
  plotData = get(paste("eKEGG",module,sep="_"))
  a=barplot(plotData, drop=TRUE, showCategory = 10, 
            title = paste0("KEGG enrichment-",module),font.size=8)+
    ggplot2::theme(legend.position =  "none",
                   plot.title = element_text(face="bold",hjust=0.5,vjust=0.5))
  print(a)
}
dev.off()

# GO graph
for (module in modulesOfInterest){
  pdf(file = paste0("GO graphs-",module,".pdf"))
  for (ont in ontologies){
    plotData = get(paste("ego",module,ont,sep="_"))
    a=plotGOgraph(plotData)
    print(a)
  }
  dev.off()
}

# Dotplots of GO enrichment
pdf(file = paste0("Dotplots of GO enrichment.pdf"),height = 8, width = 20)
for (module in modulesOfInterest){
  for (ont in ontologies){
    plotData = get(paste("ego",module,ont,sep="_"))
    a=dotplot(plotData, showCategory = 30, label_format=10) + 
      ggplot2::labs(title=ifelse(ont=="BP", paste0("Biological Pathway-",module),
                     ifelse(ont=="CC",paste0("Cellular Component-",module),
                            paste0("Molecular Function-",module))))+
      ggplot2::theme(plot.title=element_text(face="bold",hjust=0.5,vjust=0.5),
                     legend.position =  "none")
    print(a)
  }
}
dev.off()

# Dotplots of KEGG enrichment
pdf(file = paste0("Dotplots of KEGG enrichment.pdf"),height = 8, width = 20)
for (module in modulesOfInterest){
  plotData = get(paste("eKEGG",module,sep="_"))
    a=dotplot(plotData, showCategory = 30, label_format=10) + 
      ggplot2::labs(title=paste0("KEGG enrichment-", module))+
      ggplot2::theme(plot.title=element_text(face="bold",hjust=0.5,vjust=0.5),legend.position =  "none"))
    print(a)
  }
dev.off()

#CNET plot
# cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms
# or KEGG pathways) as a network (helpful to see which genes are involved in
# enriched pathways and genes that may belong to multiple annotation categories).

# convert gene ID to Symbol
for (module in modulesOfInterest){
  for (disease in diseases){
    source = get(paste("fcData",module,disease,sep="_"))
    outputVar = paste("symbol",module,disease,sep="_")
    symbols = convertToSymbol(row.names(source))
    assign(outputVar, source[source$EnsemblID %in% (symbols)[,1],2])
    dummy2 = get(outputVar)
    names(dummy2) = (symbols)[,2]
    dummy2 = sort(dummy2,decreasing = TRUE)
    assign(outputVar,dummy2)
  }
}

# cnet plot for GO
pdf(file = paste0("CNET plot of GO terms.pdf"),height = 8, width = 20)
for (module in modulesOfInterest){
  #for (ont in ontology){
    plotData = get(paste("ego",module,"ALL",sep="_"))
    y = plotData[plotData$qvalue < 0.1, asis=T] #
    for (disease in diseases){
      foldChange = get(paste("symbol",module,disease,sep="_"))
      ego1<-setReadable(y,'org.Hs.eg.db','ENSEMBL')
      a=cnetplot(ego1, foldChange=foldChange,categorySize="pvalue", cex_label_category=1, cex_label_gene = 1.1) + 
        ggraph::geom_edge_link(width=0.5, alpha = 0.1) + 
        ggplot2::labs(title=paste("Genes in enriched GO terms with qvalue<0.1",module,disease,sep="-"))+
        ggplot2::theme(plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5),
                       legend.title = element_text(size=13,face="bold"))
      print(a)
    }
  }
#}
dev.off()

# Need to figure out how to change node color. WHICH FOLD CHANGE TO USE? - Fold changes are significantly different between diseases
# CNET plot for kegg
pdf(file = paste0("CNET plot of KEGG terms.pdf"),height = 8, width = 20)
for (module in modulesOfInterest){
  plotData = get(paste("eKEGG",module,sep="_"))
  y = plotData[plotData$qvalue < 0.1, asis=T] #
  if (dim(y)[1]>0){
  for (disease in diseases){  
    foldChange = get(paste("symbol",module,disease,sep="_"))
    ego1<-setReadable(y,'org.Hs.eg.db','ENTREZID')
    a=cnetplot(ego1, foldChange=foldChange,categorySize="qvalue", cex_label_category=1, cex_label_gene = 1.1 ) + 
      ggraph::geom_edge_link(width=0.5, alpha = 0.1) + 
      ggplot2::labs(title=paste("Genes in enriched KEGG terms with qvalue<0.1",module,disease,sep="-"))+
      ggplot2::theme(plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5),
                     legend.title = element_text(size=13,face="bold"))
    print(a)
  }
  }
}
dev.off()

# categorySize can be scaled by 'pvalue' or 'geneNum'
p1 <- cnetplot(ego1, categorySize="qvalue", foldChange=symbol_darkred_DCM)
p2 <- cnetplot(ego1, foldChange=symbol_darkred_DCM, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p2, p3, ncol=2, labels=LETTERS[1:3], rel_widths=c(1.2, 1.2))

# heat plots of genes involved in the top 30 GO terms
pdf(file = paste0("Heatplot of genes involved in enriched GO terms.pdf"),height = 8, width = 20)
for (module in modulesOfInterest){
  for (disease in diseases){
    #for (ont in ontology){
      plotData = get(paste("ego",module,ont,sep="_"))
      y = plotData[plotData$qvalue < 0.1, asis=T]
      if (dim(y)[1]>0){
      foldChange = get(paste("symbol",module,disease,sep="_"))
      ego1<-setReadable(y,'org.Hs.eg.db','ENSEMBL')
      a=heatplot(ego1, foldChange=foldChange,showCategory = 30) + 
        ggplot2::labs(title=paste("Genes involved in enriched GO terms with qvalue<0.1 ",module,disease,sep="-"))+
        ggplot2::theme(plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5),
                       legend.title = element_text(size=13,face="bold"))
      print(a)
      }
    }
  }
#}
dev.off() 

# heat plots of genes involved in the top 5 KEGG terms
pdf(file = paste0("Heatplot of genes involved in enriched KEGG terms.pdf"),height = 8, width = 20)
for (module in modulesOfInterest){
  for (disease in diseases){
      plotData = get(paste("eKEGG",module,sep="_"))
      y = plotData[plotData$qvalue < 0.1, asis=T]
      if (dim(y)[1]>0){
      foldChange = get(paste("symbol",module,disease,sep="_"))
      ego1<-setReadable(plotData,'org.Hs.eg.db','ENTREZID')
      a=heatplot(ego1, foldChange=foldChange,showCategory = 10) + 
        ggplot2::labs(title=paste("Genes involved in enriched KEGG terms with qvalue<0.1",module,disease,sep="-"))+
        ggplot2::theme(plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5),
                       legend.title = element_text(size=13,face="bold"))
      print(a)
      }
    }
  }
dev.off() 

# dotplot for all three ontologies
pdf(file = "Dotplots for all GO ontologies.pdf", width = 10)
for (module in modulesOfInterest){
  source = get(paste0("ego_",module,"_ALL"))
  y = source[source$qvalue < 0.1, asis=T]
  a = dotplot(y,split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")+
    ggplot2::labs(title=paste0("Dotplot of highly enriched terms by ontology-", module))+
    ggplot2::theme(plot.title=element_text(face="bold",hjust=0.5,vjust=0.5))
  print(a)
}
dev.off()

# DAG plot - ugly plots
pdf(file = paste0("GO DAG plots.pdf"),height = 8, width = 20)
for (module in modulesOfInterest){
  for (ont in ontologies){
  plotData = get(paste("ego",module,ont,sep="_"))
  a=goplot(plotData,showCategory = 10,color = "p.adjust", layout = "sugiyama",geom = "text")+
    ggplot2::labs(title = paste("DAG plot",module,ont,sep="-")) +
    ggplot2::theme(legend.text=element_text(size=10),
                   legend.title=element_text(size=13, face = "bold"),
                   plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5))
  print(a)
  }
}
dev.off()
  
# Emphasizes the genes overlapping among different gene sets.
pdf(file = paste0("UpSet plots.pdf"),height = 14, width = 20)
for (module in modulesOfInterest){
  for (ont in ontology){
  plotData = get(paste("ego",module,ont,sep="_"))
  a=upsetplot(plotData,n=30)+
    ggplot2::labs(title = paste("UpSet plot of top 30 terms",module,ont,sep="-")) +
    ggplot2::theme(legend.text=element_text(size=10),
                   legend.title=element_text(size=13, face = "bold"),
                   plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5))
  print(a)
  }
}
dev.off()

# Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way,
# mutually overlapping gene sets tend to cluster together, making it easy to identify functional modules.

pdf(file="Enrichment map plots.pdf", height=10, width=10)
for (module in modulesOfInterest) {
  for (ont in ontologies) {
    source = get(paste("ego",module,ont,sep="_"))
    y = source[source$qvalue < 0.1, asis=T]
    a=emapplot(pairwise_termsim(y),showCategory = 60,color = "qvalue",
             layout = "nicely", min_edge = 0.5,cex_label_category = 0.5,
             cex_category = NULL,cex_line = 0)+ 
      ggplot2::labs(title = paste("Enrichment map",module,ont,sep="-"))+
      ggplot2::theme(legend.text=element_text(size=10),
                     legend.title=element_text(size=13, face = "bold"),
                     plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5))+
      ggraph::geom_edge_link(width=0.5, alpha = 0) 
    print(a)
    }
}
dev.off()

#-------------------------------------------------------------------------------#
# 9 - GSEA (GO and KEGG) of WGCNA modules - clusterprofiler
#-------------------------------------------------------------------------------

# GSEA using GO
for (module in modulesOfInterest){
  for (disease in diseases){
    file = paste("GSEresults",module,disease,sep="_",".xlsx")
    for (ont in ontologies){
    source <- get(paste("rankedGeneList",module,disease,sep="_"))
    dummy1 <- sort(source, decreasing= TRUE)
    source <- dummy1
    outputVar <- paste("GSEA",module,disease,ont,sep="_")
    assign(outputVar,gseGO(geneList = source,
                           OrgDb = org.Hs.eg.db,
                           ont = ont,
                           keyType = "ENSEMBL",
                           minGSSize = 2,
                           maxGSSize = 500,
                           pvalueCutoff = 1,
                           verbose = FALSE))
    if ((dim(get(outputVar))[1])!=0){
      write.xlsx(get(outputVar),
               file = file,
               sheet = ifelse(ont == "BP", "BiologicalProcess",
                              ifelse(ont == "CC", "CellularComponent","MolecularFunction")),
               row.names = FALSE,
               append = TRUE)
    }
    }
  }
}

# GSEA using KEGG - Only darkred_DCM, darkred_HCM, darkred_PPCM, lightyellow_DCM have significant results

for (module in modulesOfInterest){
  for (disease in diseases){
    file = paste("GSEresultsKEGG",module,disease,sep="_",".xlsx")
    source <- get(paste("entrezRanked",module,disease,sep="_"))
    outputVar <- paste("GSEAkegg",module,disease,sep="_")
    assign(outputVar,gseKEGG(source, organism = "hsa", keyType = "ncbi-geneid", 
                             exponent = 1,minGSSize = 2, maxGSSize = 500, eps = 1e-10,
                             pvalueCutoff = 0.5, pAdjustMethod = "BH", verbose = TRUE,
                             use_internal_data = FALSE, seed = FALSE, by = "fgsea"))
    if ((dim(get(outputVar))[1])!=0){
      write.xlsx(get(outputVar),
                 file = file,
                 row.names = FALSE,
                 append = TRUE)}
    }
}

# plotting ridgeplot 
pdf(file = "ridgeplots of GSEA_kegg.pdf")
for (module in modulesOfInterest){
  for (disease in diseases){
    source = get(paste("GSEAkegg",module,disease,sep="_"))
    y = source[source$qvalues < 0.1, asis=T]
    if ((dim(y)[1])!=0){
      a = ridgeplot(y, showCategory = 30, fill = "qvalues", core_enrichment = TRUE, label_format = 30)+
        ggplot2::labs(title=paste("GSEA Ridge plot",module,disease,sep="-"))+
        ggplot2::theme(plot.title=element_text(face="bold",hjust=0.5,vjust=0.5))
      print(a)
    }
  }
}
dev.off()

# visualizes the pathway and saves as png in working directory

#2099/51196/6567/4624/1734 - genes in kegg identified thyroid pathway
#ENSG00000104324/ENSG00000137857/ENSG00000140254/ENSG00000143387/ENSG00000147100/ENSG00000211448 - GO identified thyroid pathway
a=entrezRanked_tan_DCM[c("2099","51196","6567","4624","1734")]
b=convertToEntrez(c("ENSG00000104324","ENSG00000137857",
                    "ENSG00000140254","ENSG00000143387","ENSG00000147100","ENSG00000211448"))
pathview(gene.data = a, pathway.id = "hsa04919", species = "human",map.null=FALSE)


# Gsea with wikipathways - no term enriched?! - ordered list of fc and entrez id

for (module in modulesOfInterest){
  for (disease in diseases){
    file = paste("GSEresultsWP",module,disease,sep="_",".xlsx")
    entNames <- convertToEntrez(get(paste0("genesIn_",module)))
    source <- get(paste("fcData",module,disease,sep="_"))
    column <- paste0(disease,"fc")
    dummy1<- source[source[,1] %in% entNames[,1],2]
    names(dummy1) <- entNames[,2]
    source <- sort((dummy1),decreasing=TRUE)
    outputVar <- paste("GSEAwp",module,disease,sep="_")
    assign(outputVar,gseWP(source, "Homo sapiens",
                           pvalueCutoff = 0.5,
                           pAdjustMethod = "BH",
                           #qvalueCutoff = 0.5,
                           minGSSize = 2,
                           maxGSSize = 500))
      if ((dim(get(outputVar))[1])!=0){
      write.xlsx(get(outputVar),
                 file = file,
                 row.names = FALSE,
                 append = TRUE)}
    
  }
}

# Enriched KEGG Modules - characteristic gene sets that can be linked to 
# specific metabolic capacities and other phenotypic features, so that they can
# be used for automatic interpretation of genome and metagenome data
# Still have to figure out why this doesn't give any results.
#mkk <- enrichMKEGG(gene = source, organism = "hsa", keyType = "ncbi-geneid")
for (module in modulesOfInterest){
  source = get(paste0("genesIn_",module))
  dummy1 = convertToEntrez(source)
  source = dummy1[,2]
  file = paste0("MKEGG_EnrichmentResults_",module,".xlsx")
  outputVar = paste("eMKEGG",module,sep="_")
  assign(outputVar, enrichMKEGG(gene = source, organism = "hsa", keyType = "ncbi-geneid",
                               pvalueCutoff = 0.5,
                               pAdjustMethod = "BH",
                               universe = EntrezIDs[,2],
                               minGSSize = 2,
                               maxGSSize = 500,
                               qvalueCutoff = 0.5))
  resKEGG <- get(outputVar)
  if (!is.null(resKEGG)){
  write.xlsx(resKEGG,
             file = file,
             row.names = FALSE)}
}


# gsea with kegg modules. 

for (module in modulesOfInterest){
  if (module != "darkred"){ #Darkred module does not work for some reason - no ids can be mapped.
  for (disease in diseases){
    file = paste("GSEresultsMKEGG",module,disease,sep="_",".xlsx")
    entNames <- convertToEntrez(get(paste0("genesIn_",module)))
    source <- get(paste("fcData",module,disease,sep="_"))
    column <- paste0(disease,"fc")
    dummy1<- source[source[,1] %in% entNames[,1],2]
    names(dummy1) <- entNames[,2]
    source <- sort((dummy1),decreasing=TRUE)
    outputVar <- paste("GSEAmkegg",module,disease,sep="_")
    assign(outputVar,gseMKEGG(geneList = source, 
                              keyType = "ncbi-geneid", 
                              organism = "hsa",
                              pvalueCutoff = 0.5,
                              pAdjustMethod = "BH",
                              minGSSize = 2,
                              maxGSSize = 500))
   if ((dim(get(outputVar))[1])!=0){
   write.xlsx(get(outputVar),
             file = file,
             row.names = FALSE,
             append = TRUE)}
  }
    
  }
}

#-------------------------------------------------------------------------------#
# 10 - Biological theme comparison of a few modules
#-------------------------------------------------------------------------------
#compareCluster() requires named list of gene ids of different modules. 
#each element of the list must be vector of gene ids. IDs need to be Entrez
# fun = "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPath-way"

# list containing Entrez ids of top three modules. 
entList <- list()
for (module in modulesOfInterest){
  geneInputList <- get(paste0("genesIn_",module))
  a = convertToEntrez(geneInputList)[,2]
  entList[[paste0(module)]] = a
}

# list containing Ensembl ids of top three modules. 
ensList <- list()
for (module in modulesOfInterest){
  ensList[[paste0(module)]] <- get(paste0("genesIn_",module))
}

# comparing the enriched terms between clusters
clusterComparison=compareCluster(geneCluster = entList,fun="enrichGO",
                                 OrgDb = org.Hs.eg.db,pvalueCutoff=0.05, 
                                 pAdjustMethod="BH",ont="BP",readable=T)
b=data.frame(clusterComparison)
dotplot(clusterComparison,includeAll=F)

#cnetplot(setReadable(clusterComparison,OrgDb = org.Hs.eg.db), categorySize="pvalue", cex_label_category=1, cex_label_gene = 1.1 ) + 
#  ggraph::geom_edge_link(width=0.5, alpha = 0.1) + 
#  ggplot2::labs(title="Genes ")+
#  ggplot2::theme(plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5),
#                 legend.title = element_text(size=13,face="bold"))

a=emapplot(pairwise_termsim(clusterComparison),showCategory = 60,color = "qvalue",
           layout = "nicely", min_edge = 0.5,cex_label_category = 0.8,
           cex_category = 0.8,cex_line = 0)+ 
  ggplot2::labs(title = "Enrichment map for comparison of clusters",
                fill = "Module")+
  ggplot2::theme(legend.text=element_text(size=10),
                 legend.title=element_text(size=13, face = "bold"),
                 plot.title = element_text(size=20, face = "bold",hjust=0.5,vjust=0.5))+
  ggraph::geom_edge_link(width=0.5, alpha = 0) 
print(a)



# Don't really know how to interpret this. It allows us to create one or more groups based on 
# FC values of genes within a module. 
mydf <- data.frame(Entrez=names(entrezRanked_darkred_HCM), FC=entrezRanked_darkred_HCM)
mydf <- mydf[abs(mydf$FC) > log2(1.5),]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"
formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")
head(as.data.frame(formula_res))
dotplot(formula_res)


#-------------------------------------------------------------------------------#
# 11 - topGO on modules identified by WGCNA 
#-------------------------------------------------------------------------------

for (module in modulesOfInterest){
  InterestingGenes <- get(paste0("genesIn_",module))
  for (ont in ontologies) {
    outputVar <- paste("allRes",module,ont,sep="_") #setting output variable name
    file <- paste("TopGO_Results_",module,".xlsx", sep="") #setting file name
    print(ont)
    dataForTopGO <- annFUN.org(ont, mapping = "org.Hs.eg.db", ID = "ensembl")
    universeGenes <- unique(unlist(dataForTopGO))
    factorGenesTopGO <- factor(as.integer(universeGenes %in% InterestingGenes))
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


#-------------------------------------------------------------------------------#
# 12 - PathfindR on modules identified by WGCNA
#-------------------------------------------------------------------------------

#running pathfindR on all modules - PIN=STRING, GSET=KEGG. 
for (module in c("darkred","tan")){# modulesOfInterest){
  for (disease in diseases){
    dirEnd <- paste("pathfindr",module,disease,sep="_")
    dir <- paste0("C:/Users/kaila/Documents/MSc/Period 2/Experimental methods and data management/R skill sessions/Data/pathfindR_PINstring_genesetKEGG/",dirEnd)
    source <- get(paste("fcData",module,disease,sep="_"))
    dummy1 <- source[,1]
    dummy1 <- convertToSymbol(dummy1)
    dummy2 <- cbind(dummy1,source[source$EnsemblID %in% dummy1[,1],])
    source <- dummy2[,c(2,4,5)]
    outputVar <- paste("pathRes",module,disease,sep="_")
    assign(outputVar, run_pathfindR(source,
                                    p_val_threshold = 0.05,
                                    convert2alias = TRUE,
                                    pin_name_path = "STRING", # Name of the chosen PIN - "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG"
                                    sig_gene_thr = 0.02, # Threshold for minimum proportion of significant genes used in filtering active subnetworks
                                    gene_sets = "KEGG", # gene set used for enrichment analysis - "KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC", "GO-MF"
                                    min_gset_size = 10,
                                    max_gset_size = 500,
                                    output_dir = dir,
                                    plot_enrichment_chart = TRUE))
  }
}

#performing hierarchical clustering and saving plots
pdf(file = "Pathfinder Cluster enrichment plots.pdf", height=9, width = 9)
for (module in modulesOfInterest){
  for (disease in diseases){
    source <- get(paste("pathRes",module,disease,sep="_"))
    outputVar <- paste("clustered",module,disease,sep="_")
    chart <- paste("enrichmentplotByCluster",module,disease,sep="_")
    if (dim(source)[1]!=0){
      assign(outputVar,cluster_enriched_terms(source,plot_dend = FALSE,
                                              plot_clusters_graph = FALSE))
      
      assign(chart,enrichment_chart(get(outputVar), plot_by_cluster = TRUE)+
               ggplot2::labs(title = paste("Cluster enrichment plot",module,disease,sep="-"))+
               ggplot2::theme(legend.text=element_text(size=10),
                              legend.title=element_text(size=13, face = "bold"),
                              plot.title = element_text(size=17, face = "bold",hjust=0.5,vjust=0.5)))
      print(get(chart))
    }
  }
}
dev.off()

# term-gene graphs for pathfindr term-clusters 
pdf(file = "Pathfinder Cluster term-gene plots.pdf", height=10, width = 10)
for (module in modulesOfInterest){
  for (disease in diseases){
    #if (module != "grey60"){ #grey60 has no clustered terms so it throws an error. 
    source <- get(paste("clustered",module,disease,sep="_"))
    chart2 <- paste("ClusterTGpathfindr",module,disease,sep="_")
    if (dim(source)[1]!=0){
     assign(chart2,term_gene_graph(source,node_size = "num_genes",use_description = TRUE, layout = "nicely")+
               ggplot2::labs(title = paste("Cluster term-gene plot",module,disease,sep="-"), fill = "Number of genes")+
               ggplot2::theme(legend.text=element_text(size=10),
                              legend.title=element_text(size=13, face = "bold"),
                              plot.title = element_text(size=17, face = "bold",hjust=0.5,vjust=0.5)))
      print(get(chart2))
    }
    }
  }
#}
dev.off()

# creating input for next step - plotting heatmaps of enriched terms vs genes involved. 
for (module in modulesOfInterest){
  for (disease in diseases){
    outputVar <- paste("forVisualizePathf",module,disease,sep="_")
    source <- get(paste("fcData",module,disease,sep="_"))
    dummy1 <- source[,1]
    dummy1 <- convertToSymbol(dummy1)
    dummy2 <- cbind(dummy1,source[source$EnsemblID %in% dummy1[,1],])
    assign(outputVar,dummy2[,c(2,4,5)])
  }
}

#plotting heatmap of enriched terms vs genes involved - FAILS FOR TAN AND LIGHTYELLOW
pdf(file="heatmap - enriched terms vs genes involved.pdf", height = 10, width = 10)
for (module in modulesOfInterest){
  for (disease in diseases){
    source <- get(paste("pathRes",module,disease,sep="_"))
    genesInput <- get(paste("forVisualizePathf",module,disease,sep="_"))
    chart <- paste("term-geneHeatmap",module,disease,sep="_")
    if(dim(source)[1]!=0){
    assign(chart,
           term_gene_heatmap(result_df= source, use_description = TRUE,
                             high = "green", mid = "black",low = "red")+
             ggplot2::labs(title = paste("Enriched terms vs genes involved",module,disease,sep="-"), fill = "Fold Change")+
             ggplot2::theme(legend.text=element_text(size=10),
                            legend.title=element_text(size=13, face = "bold"),
                            plot.title = element_text(size=17, face = "bold",hjust=0.5,vjust=0.5)))
    print(get(chart))
    }
  }
}
dev.off()


# plot term-gene graphs
for (module in modulesOfInterest){
  for (disease in diseases){
    source <- get(paste("pathRes",module,disease,sep="_"))
    if(dim(source)[1]!=0){
    chart2 <- paste("term-geneGraph",module,disease,sep="_")
    assign(chart2,term_gene_graph(source,node_size = "p_val"))
    print(get(chart2))
    }
  }
}

# upset plots
pdf(file="Upset plots_pathfindrResults.pdf", height = 14, width =14)
for (module in modulesOfInterest){
  for (disease in diseases){
    source <- get(paste("pathRes",module,disease,sep="_"))
    genesInput <- get(paste("forVisualizePathf",module,disease,sep="_"))
    chart3 <- paste("upset",module,disease,sep="_")
    if(dim(source)[1]!=0){
    assign(chart3, UpSet_plot(source,num_terms = 10,use_description = TRUE,genes_df = genesInput,high = "green", mid = "black",low = "red"))
    print(get(chart3))
    }
  }
}
dev.off()

# Vector of "Case" IDs
casesPPCM <- row.names(sampleInfo[sampleInfo$Disease == "PPCM", ])
casesDCM <- row.names(sampleInfo[sampleInfo$Disease == "DCM", ])
casesHCM <- row.names(sampleInfo[sampleInfo$Disease == "HCM", ])

# Calculate scores for representative terms and plot heat map using term descriptions - 
#This allows the user to individually examine the scores and infer how a term 
#is overall altered (activated or repressed) in a given sample or a group of samples
pdf(file="Heatmap of score matrix of enriched terms per sample.pdf",height=6,width=17)
for (module in c("tan", "darkred", "lightyellow")){#modulesOfInterest){
  for (disease in diseases){
    if (!((module=="grey60"|module=="royalblue") & disease=="PPCM")){
    enrichment <- get(paste("clustered", module, disease, sep="_"))
    cases <- get(paste0("cases",toupper(disease)))
    outputVar <- paste("score_matrix",module,disease,sep="_")
    assign(outputVar,score_terms(enrichment_table = enrichment,#[enrichment$Status == "Representative",],
                                 exp_mat = as.matrix(ABgxSymbols),
                                 cases = cases,
                                 control_title = "Healthy", # default = "Control"
                                 use_description = (!(module=="tan" & disease=="PPCM")), # default FALSE
                                 label_samples = FALSE, # default = TRUE
                                 case_title = as.character(disease),# default = "Case"
                                 low = "#f7797d", # default = "green"
                                 mid = "#fffde4", # default = "black"
                                 high = "#1f4037")) # default = "red"
    }
  }
}
dev.off()

#-------------------------------------------------------------------------------#
# 13 - Compare any two module pathfindR results
#-------------------------------------------------------------------------------
combined_df <- combine_pathfindR_results(result_A = pathRes_tan_PPCM, 
                                         result_B = pathRes_darkred_PPCM, 
                                         plot_common = F)
combined_results_graph(combined_df,
                       selected_terms = "common", 
                       use_description = TRUE,
                       layout = "stress",
                       node_size = "p_val")

# 12 - Dataframe of differentially expressed genes in WGCNA modules


