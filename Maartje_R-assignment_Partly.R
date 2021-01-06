#=============================================================================#
# THIS FILE CONTAINS PART OF MY R-SCRIPT FOR MBS1005      
# 
# I excluded chunks of code which I don't think are that useful for us. 
# (e.g. I removed the biggest part of the PCA plots, because I didn't have enough 
# time in the end to make a loop, so it was just tediously long). 
# 
# Actually,  the only thing that might have some added value is the code for the 
# volcano plots in part 3
#
# (Also: I have no idea if what I did for part 5 is even remotely right.) 
#
#=============================================================================#

#=============================================================================#
# Install/Load Packages
#=============================================================================#


#----------------------------------------------------------------------------------------------------------------------------#
# 
# Note: 'Un-comment' following code if packages are not installed (and change 
#       'setRepositories()' from FALSE to TRUE if necessary).
#
# setRepositories(FALSE, 1:10)                                                                                               
# install.packages(c('ggplot2', 'limma', 'pcaMethods', 'biomaRt','vioplot', 'gtools', 'arsenal', 'tidyr', 'qvalue','dplyr')) 
#
#----------------------------------------------------------------------------------------------------------------------------#

require(ggplot2)
require(limma)
require(pcaMethods)
require(gtools)
require(vioplot)
require(biomaRt)
require(arsenal) 
require(tidyr)
require(qvalue)
require(dplyr)
options(stringsAsFactors = F)

#=============================================================================#
# Part 1: Data import
#=============================================================================#

# Setting the working directory

setwd("C:\\Users\\maart\\OneDrive\\Academic\\P2_Experimental_Design_and_Data_Management\\skills_R")


#-----------------------------------------------------------------------------#
# Part 1.1 Import the data
#-----------------------------------------------------------------------------#

# Importing the sample characteristics data

sampleInfo <- read.delim("MAGNET_SampleData_19112020.txt", as.is = T, 
                         row.names = 1)

# Importing the gene expression data in counts per million (CPM)

gxData <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", as.is = T,
                     row.names = 1)

#-----------------------------------------------------------------------------#
# Part 1.2 Export a table and/or figure(s) of sample characteristics
#-----------------------------------------------------------------------------#

## REMOVED ##

#-----------------------------------------------------------------------------#
# Part 1.3 Transform CPM values to FPKM values 
#-----------------------------------------------------------------------------#

## REMOVED ##


#=============================================================================#
# Part 2:  Diagnostic plots
#=============================================================================#

#-----------------------------------------------------------------------------#
# Part 2.1 Figures containing boxplots
#-----------------------------------------------------------------------------#
# Creating four seperate PNG files exported to the working directory 
# each file contains a boxplot for a different disease group.
# The width and height can be adapted to make the figure more readable.

png(filename = "Boxplot_Donor.png", width = 2000, height = 1000) 
boxplot(gxData[, sampleInfo$Disease == "Donor"],
        las = 2,
        xlab = "SampleID", 
        ylab = "Expression", 
        main = "Donor samples")
dev.off() 

png(filename = "Boxplots_DCM.png", width = 2000, height = 1000) 
boxplot(gxData[, sampleInfo$Disease == "DCM"],
        las = 2,
        xlab = "SampleID", 
        ylab = "Expression",
        main = "DCM samples")
dev.off()

png(filename = "Boxplot_HCM.png", width = 1500, height = 500) 
boxplot(gxData[, sampleInfo$Disease == "HCM"],
        las = 2,
        xlab = "SampleID", 
        ylab = "Expression", 
        main = "HCM samples")
dev.off()

png(filename = "Boxplots_PPCM.png", width = 500, height = 500) 
boxplot(gxData[, sampleInfo$Disease == "PPCM"],
        xlab = "SampleID", 
        ylab = "Expression", 
        main = "PPCM samples")
dev.off()


#-----------------------------------------------------------------------------#
# Part 2.2 Principle component analysis (PCA)
#-----------------------------------------------------------------------------#
# Use the function 'pca' from the 'pcaMethods' package to run a PCA 
# for 5 principal components (PCs) and save it as pca_data. 
# The gxData is transposed here to get the analysis to run over the samples.

pca_data <- pca(t(gxData), nPcs = 5)

#Now you have PCA results for all samples, which you can see by printing the 
#pca_data@scores

print(pca_data@scores)

# Additionally you can take a look at the proportion of variance (R2) 
# that is explained by the 5 PCs:
#
# The reason to inspect the R2 is because it can give some indication whether 
# the number of PCs makes sense.For example, if the R2 explained by last PC is small, 
# that might be a reason to not include in the PCA plots. 
# On the other hand, if the last PC still adds a lot to the proportion of variance explained,
# then it might be good to expand the PCA 
# to see if there are even more PCs which explain a lot of variance. 
#
# NOTE: This is just an eyeballing method!
#
# Inspect the R2 for each PC by printing an plotting: 

print(pca_data)

# Combine the pca_data with the sample characteristic information and save as plot_PCA_char

plot_PCA_char <- cbind(data.frame(pca_data@scores), sampleInfo)

# The following part of code uses this plot_PCA_characteristics data 
# to make PCA plots.For all possible PC duos three plots are made, 
# one for each covariate (Age, Sex, Ethnicity). 
#
# The plots are colored by the covariate.
# The variable of interest (Disease) is represented by different shapes. 


#---------------------------#
# PCA plots                 #
#---------------------------# 

## REMOVED ##

#---------------------------#
# Plots for PC1 against PC2 #
#---------------------------# 

# Plot colored by Age
png(filename = "PCA1vs2_age", width = 1000, height = 1000) 
ggplot(plot_PCA_char, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = Disease, col = Age))
dev.off()


# Plot colored by Sex  
ggplot(plot_PCA_char, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = Disease, col = Sex))

# Plot colored by Ethnicity
ggplot(plot_PCA_char, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = Disease, col = Ethnicity))


## etc... REMOVED ##


#=============================================================================#
# Part 3: Statistical analysis
#=============================================================================#
# The 'limma' package is used for the statistical analysis.

# create a design matrix

design <- model.matrix(~0 + Disease + Age + Sex + Ethnicity, data = sampleInfo) 

rownames(design) <- rownames(sampleInfo)


# Because there are more than 2 groups for the independent variable 'Disease'
# use 'makeContrasts' to specify the contrasts you want to test).

contrast_matrix <- makeContrasts(
  DCMvsDonor = DiseaseDCM-DiseaseDonor, 
  HCMvsDonor = DiseaseHCM-DiseaseDonor, 
  PPCMvsDonor = DiseasePPCM-DiseaseDonor,
  levels=design
)

# Fit the coefficients of the linear model with 'lmFit' by passing it the
# gene expression data and the design matrix ('gxdata' and 'design', 
# respectively).

fit <- lmFit(gxData, design)

# Compute the estimated coefficients and standard error for the contrasts included 
# in the linear model
fit_contrast <- contrasts.fit(fit, contrast_matrix)

# Use 'eBayes()' to compute the moderated t-statistics, 
# F-statistics and logg-odds of differential expression based on the empirical
# Bayes moderation of the standard errors towards a global value.

fit_contrast <- eBayes(fit_contrast, trend=TRUE)

# Generate volcano plots for each comparison to visualize differential expression
# for each pairwise contrast

for (i in 1:ncol(fit_contrast)) {
  png(
    filename = paste("plot_volcano_disease", (i), "_vs_donor.png", sep = ""), 
    width = 1000, 
    height = 1000
  )
  volcanoplot(
    fit = fit_contrast[,i], 
    style = "p-value",
    main = colnames(fit_contrast)[i], 
    col=ifelse(fit_contrast[,i]$p.value > 0.05,"red","green")
  )
  abline(v=log2(2))  
  abline(v=-log2(2)) 
  abline(a = -log10(0.05), b = 0) 
  dev.off()
} 

# Use the 'topTable()' function to generate summary statistics of the linear model. 
# The adj.P.val column in the output contains p-values adjusted for multiple testing.
# The adjusted p-values are computed with teh Benjamini and Hochberg's (BH) method 
# which controls for the false discovery rate (i.e. the proportion of false positives 
# when a test is called significant based on the p-value).
# NOTE: these adjusted p-values are otherwise known as q-values. 

modelstat_summary <- topTable(fit_contrast, number = nrow(gxData), adjust.method = "BH")

# Next create a new expression data frame that only includes the genes with a 
# differential expression level which are significant according to the adjusted
# p-values.

gxTop <- gxData[modelstat_summary$adj.P.Val < 0.05, ]


# Note: the following 2 steps are not vital to the script. However, it 
#       might come in handy if you want to combine the two data sets later.

# Adjust the summary statistics table such that its rows match the gxData

modelstat_summary <- gene_lmtable[match(rownames(gxData),rownames(gene_lmtable)),]

# do a quick sanity check to see if the gene IDs match the gxData

all(rownames(modelstat_summary) == rownames(gxData)) #TRUE

#=============================================================================#
# Part 4: Additional gene annotation 
#=============================================================================#
# In this part the 'biomaRt' package is used to add annotations for gene symbols, 
# gene names and gene descriptions to the gene expression data using a BioMart 
# database. 

# Note: If the exact titles of the desired database and dataset are known 
#       beforehand, you can connect to set database and dataset in one step,
#       with the command 'useMart()'. However, because of the fact that this 
#       only works if the titles are exactly right (thus up-to-date), I choose 
#       to connect to dataset step by step in the code below. With this approach
#       you can easily check the database and dataset titles, and it is easier 
#       to adapt if you want to use another BioMart database. 

# The annotion is going to be based on the gene identifiers in the 
# Ensembl database.The first step is to use 'listMarts()' which provides a list 
# of available databases. Here you can check the exact title that you need to 
# connect to the Enseble database. Because this list contains only a few items
# you can inspect it in your console. 

listMarts() 

# After checking the title of the database ('ENSEMBL_MART_ENSEMBL' at the time 
# of writing), the 'useMart()' function can be used to connect to the Ensembl 
# database. 

ensembl <- useMart("ENSEMBL_MART_ENSEMBL") 

# Because the Ensemble database contains multiple datasets (Each included 
# species has its own dataset) the next step is to create an object that can be
# checked for the title of the dataset for human genes. 

all_datasets <- listDatasets(ensembl)

# Update to connect to the specific human dataset

ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl) 

# Once connected to the dataset the 'getBM()' command can be used to retrieve 
# attributes of interest from the Enseml database. In addition to the connected 
# Ensembl data, this command requires three main input arguments: attributes, 
# filters and values for the filters. In order to properly use 'getBM()' we need
# to learn about the available arguments.

# Create object to inspect the available filters

all_filter <- listFilters(ensembl)

# Create object to inspect the available attributes 

all_attributes <- listAttributes(ensembl)

# Now the names of the desired attributes and filters, as retrieved from set 
# objects, can be implemented to run 'getBM()' with the Emsembl data. 
#
# The filter 'esembl_gene_id' is used to restrict the output to gene stable IDs, 
# corresponding to the IDs in the 'gxData'. This is also included as an attribute
# such that the expression data can be compared with or combined to the annotion
# information based on these IDs.
# This is followed by the attributes to retrieve gene symbols, names and descriptions, 
# respectively. The attributes provides information on the chromosome name, this
# will be used in part 5 of the script. 


annot_ensembl <- getBM(
  attributes= c("ensembl_gene_id", "hgnc_symbol", "external_gene_name","description","chromosome_name"),
  filters = "ensembl_gene_id",
  values = row.names(gxTop),
  mart = ensembl
)


# use 'intersect()' to find out which gene IDs are both in the 'gxData' and in 
# the 'annot_emsebl'.

overlap <- intersect(row.names(gxData), annot_ensembl$ensembl_gene_id)

# Use this information to create a CPM data frame only including the genes 
# of which annotation information could be obtained from Ensembl. 

gxTop_select <- data.frame(gxTop[overlap, ])

#=============================================================================#
# Part 5: Relative expression level of genes
#=============================================================================#
# This part of the script is used to assess which genes are expressed above
# background (i.e. noise) level. This is accomplished by comparing the 
# gene expression at the y-chromosome in females to the other gene expression in 
# females. Any y-chromosomal gene expression in female samples is noise. 
# Hence,only genes exceeding this level of expression are considered to be expressed 
# above background level. 

# First, transform the CPM values to FPKM values using the function, defined earlier.
# to control for the difference in length (because we are comparing expression 
# between different genes.   

fpkm_Top <- cpm2fpkm(gxTop_select)

# filter the genetic data to only include expression data from the female samples

fpkm_F_Top <- fpkm_Top[, sampleInfo$Sex == "Female"]


# Identify and separate the y-chromosome genes from the other genes, 
# and calculate the mean expression for the y chromosome genes use as the 
# background level to compare the other means with.


nonY_fpkmset <- fpkm_F_Top[annot_ensembl$chromosome_name != "Y", ]

Y_fpkmset <- fpkm_F_Top[annot_ensembl$chromosome_name == "Y", ]

background <- sum(rowMeans(Y_fpkmset[,]))/nrow(Y_fpkmset)

# Find the mean fpkm for the other genes, compare to the background and add to 
# annotation information

nonY_means <- data.frame(rowMeans(nonY_fpkmset[,]))

# make a new data frame for the cpm expression values, excluding all the data 
#which does not exceed background according to the fpkm Y comparison. 

exceed_background <- gxTop_select[nonY_means > background, ]

#=============================================================================#
# Part 6: Export final results
#=============================================================================#

## REMOVED ##

