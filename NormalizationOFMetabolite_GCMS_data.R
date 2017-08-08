########## Normalization of High Dimensional Metabolite MS Data ##########
################################################################################
# Normalization Liver Metabolites.
# Modified code from Dan Gatti.
# Anji Trujillo
# etrujillo2@wisc.edu
# July 29, 2017
################################################################################

##############################
# Load and install packages. #
##############################
source("https://bioconductor.org/biocLite.R") # opens up link to download bioconductor packages

biocLite("sva") # download SVA package
browseVignettes("sva")
biocLite("pcaMethods") # download pcaMethods
biocLite("limma") # supportive package for pcaMethods
biocLite("genefilter") # supportive package for comBat
biocLite(c("devtools","Biobase","bladderbatch","snpStats")) # supportive packages

library(tidyverse)
library(pcaMethods)
library(sva)
library(pamr)
library(limma)
library(genefilter)
library(devtools)
library(Biobase)
library(bladderbatch)
library(snpStats)

#####################
# Working Directory #
#####################
# Set input and output directories

getwd()
setwd("C:/Users/etrujillo/Desktop/DOProjectFolder")
input.dir = "C:/Users/etrujillo/Desktop/DOProjectFolder"
output.dir = "C:/Users/etrujillo/Desktop/DOProjectFolder/NormalizedThroughSAVandComBat"

#####################
# Load in the data. #
#####################

options(stringsAsFActors = F)
DOLiver01Aug2017 <- read.csv("01Aug2017DOPlasmaMetabolites_RAW_NoTransformation_Tier4_ForR.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE) # raw untransformed DO Liver

############################################
# Restructure data to include annotations. #
############################################

annot <- read.csv("DOWave1through4_Covariates.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE) # Read in the sample annotation.
rownames(annot) <- annot$Mouse.ID
metab = right_join(annot, DOLiver01Aug2017, by = "Mouse.ID") # Merge the sample annotation and data by mouse id.
x <- metab[grepl("Redo", metab$Mouse.ID),] # 14 samples with REDO in identifier
metab <- metab[!grepl("Redo", metab$Mouse.ID),] # removing all the rows with Redo in identifier
row.names(metab) <- metab[,1] # set row names
dim(metab) # 384 x 385

data  = as.matrix(metab[,-(1:22)]) # Split up the sample annotation from the data and convert the data into a numeric matrix.
data[data == 1] <- NA # replace all 1 with NA
#ctrl = which(annot$Mouse.ID == "Control") # determine control samples
#data = data[-ctrl,] # Remove control samples from data matrix
#annot = annot[-ctrl,] # Remove control samples from annot data
dim(data) # 384 mice X 363 features

##################
# Plotting data. #
##################
# pcaMethods wants samples in rows and variables in columns.

palette(c("mediumorchid2","mediumturquoise","olivedrab3", "darkgoldenrod1", 
          "hotpink3", "red2", "steelblue2", "sienna2","slategray4", 
          "deepskyblue", "orangered", "midnightblue"))

data.log2 = log2(data)
pc.data = pca(data.log2, method = "bpca", nPcs = 20) # iterative Bayesian model that handels missing values "NA"
pc.data.RAW = pca(data, method = "bpca", nPcs = 20) # should PCA be done with log2 transformed data or RAW??
str(pc.data)
y <-completeObs(pc.data)
z <- completeObs(pc.data.RAW)

pdf("Figures/Liver_metabolites_Tier4_Log2RAW_DOMice_PCA_20170808.pdf") # save pdf in Figures folder
  
  batch.colors = as.numeric(factor(metab$Batch.x))
  plot(scores(pc.data), pch = 16, col = batch.colors, main = "Un-normalized Liver Metabolites, Colored by Batch 20170802")
  text(scores(pc.data)[,1], scores(pc.data)[,2], labels = rownames(data.log2), 
       col = batch.colors)
dev.off()

prop.missing = rowMeans(is.na(data)) # mean of each row where there is missing data
sum(prop.missing > 0.25) # determine what samples have more than 25% of missing data # there are two samples
rownames(data)[prop.missing > 0.25] # [1] "DO109" "DO134"
keep = which(prop.missing < 0.25) # samples to be kept; all data that has less than 25% of missing data
data = data[keep,] # subset data to keep
metab = metab[keep,]
annot = as.data.frame(metab[,1:6])


str(data) #Raw data were samples contain has more than 75% of data present. two dimnames 382 samples with 363 metabolites
str(annot) #382 with 6 variables
str(metab) #382 obs with 385 variables

pdf("Figures/VariationbyLoading_LiverMetabolites_20170808.pdf", useDingbats = FALSE) # save pdf in Figures folder

  pc.data = pca(data.log2, method = "bpca", nPcs = 20)
  plot(pc.data)
  abline(h = 0.95, col = 2)
  
dev.off()

pdf("Figures/LoadingUnnomarlizedLiverMetabolites20170808.pdf", useDingbats = FALSE) # save pdf in Figures folder
  plot(loadings(pc.data), pch = 19, col = "mediumorchid2", main = "Loadings Un-normalized Liver Metabolites, 363 Metabolites") #plot the loadings
  text(loadings(pc.data)[,1], loadings(pc.data)[,2], labels = colnames(data.log2), 
       col = "midnightblue")
dev.off()
dim(loadings(pc.data)) #363 x 20

######################################################
# PCA plots of Unnormalized data by batch, sex, etc. #
######################################################

sex = factor(metab$sex) # create a factor for sex

pdf("figures/liver_metabolites_unnormalized_PCA_20170808.pdf", useDingbats = FALSE) # compiles the following plot into one pdf

  plot(scores(pc.data), pch = 16, col = as.numeric(sex),
       main = "Un-normalized Liver Metabolites Colored by Sex")
  legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
  
  batch = factor(metab$Batch.x)
  plot(scores(pc.data), pch = 16, col = as.numeric(batch),
       main = "Un-normalized Liver Metabolites Colored by Batch")
  legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
         y.intersp = 0.7)
  
  wave = factor(metab$Wave)
  plot(scores(pc.data), pch = 16, col = as.numeric(wave),
      main = "Un-normalized Liver Metabolites Colored by Wave")
  legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)))
  
  #diet.days = factor(annot$diet_days, levels = sort(unique(annot$diet_days)))
  #diet.colors = rainbow(length(levels(diet.days)) - 1)
  #plot(scores(pc.data), pch = 19, col = diet.colors[diet.days],
       #main = "Un-normalized Liver Metabolites Colored by Diet Days")
  
dev.off()

######################################
# Set up batch and model for comBat. #
######################################
# ComBat wants the data with variable in rows and samples in columns.
# use prior.plots = TRUE to give prior plots with kernel estimate of the empirical batch effect as well 

mod = model.matrix(~sex, data = metab) #num arry with 382 samples
batch = metab$Batch.x # 382


chg = 1e6
iter = 1
repeat( {

  print(paste("Iteration", iter))

  # Impute missing data.
  miss = which(is.na(data))
  print(paste(length(miss), "missing points."))
  pc.data = pca(data, method = "bpca", nPcs = 7)
  data.compl = completeObs(pc.data)

  # Batch adjust.
  pdf("figures/ComBatSignificantAnalysis_20170804.pdf", useDingbats = FALSE) # compiles the following plot into one pdf)
  data.cb = ComBat(dat = t(data.compl), batch = batch, mod = mod, prior.plots = TRUE) # ComBat normalization function
  dev.off()
  
  data.cb = t(data.cb)

  # Calculate the change.
  chg = sum((data.compl[miss] - data.cb[miss])^2)
  print(paste("   SS Change:", chg))

  # Put the missing data back in an impute again.
  if(chg > 1 & iter < 20) # if chg is greater than 1 and iter is less than 20
    { 
      data.cb[miss] = NA
      data.comBat = data.cb
      iter = iter + 1    
    }else 
    {
      data.comBat = data.cb # assign data.cb(comBat analysis) to data.log
      break
    }
  }
)

########################################
# Remove or average duplicate samples. #
########################################

dupl = which(duplicated(rownames(data.comBat)))  # no names are duplicated
dupl.data = data[rownames(data.comBat) %in% rownames(data.comBat)[dupl],] # index rows where there are duplicates
stopifnot(rownames(data) == rownames(data.comBat))
prop.missing = rowMeans(is.na(data))
unique.samples = unique(rownames(data.comBat)) #382 unique samples
keep = rep(FALSE, nrow(data.comBat))

for(i in 1:length(unique.samples)) 
  {
    sample = unique.samples[i]
    wh = which(rownames(data.comBat) == sample)
    wh = wh[which.min(prop.missing[wh])]
    keep[wh] = TRUE
  }

data.comBat.keep = data.comBat[keep,]
metab.match = metab[match(rownames(data.comBat.keep), metab$Mouse.ID),] #return the position of first occerences of the vector1 in vector2
annot.match = annot[match(rownames(data.comBat.keep), annot$Mouse.ID),]
x <- match(rownames(data.comBat.keep), annot$Mouse.ID)
y <- match(rownames(data.comBat.keep), metab$Mouse.ID)

# Merge in the Chr M and Y info.
# attie_MY = read_csv(paste0(input.dir, "attie_sample_info_ChrM_Y.csv"))
# annot = right_join(annot, attie_MY, by = "Mouse.ID")
# annot = annot[,c(1:10, 13:15)]
# colnames(annot) = sub("\\.x", "", colnames(annot))

data.comBat.withcorrectedMouseID = data.frame(Mouse.ID = rownames(data.comBat.keep), data.comBat.keep) 
data.out = right_join(annot, data.comBat.withcorrectedMouseID, by = "Mouse.ID") 
row.names(data.out) <- data.out$Mouse.ID
str(data.out) # data frame with 382 obst and 369 variables
saveRDS(data.out, file = paste0(output.dir, "EAT_Liver_GCMS_Metabolites_ComBatNormalized.rds"))

#################
# Z-score data. # 
#################

# Transform each analyte into Z-scores.
data.rz = data.out

rankZ = function(x) 
  {
    x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
    return(qnorm(x))
  } # rankZ()

 for(i in 14:ncol(data.rz)) 
  {
    data.rz[,i] = rankZ(data.rz[,i])
  }

# saveRDS(data.rz, file = paste0(output.dir, "attie_liver_metabolites_zscore_normalized.rds"))

####################################################
# PCA plots of Normalized data by batch, sex, etc. #
####################################################
 
x <- data.out[,-(1:6)]
 
pdf("figures/liver_metabolites_normalized_PCA_20170807.pdf", width = 12, height = 7, useDingbats = FALSE)

  pc.data = pca(as.matrix(data.out[,-(1:6)]), method = "bpca", nPcs = 20)
  
  layout(matrix(1:2, 1, 2))
  
  sex = factor(data.out$sex)
  plot(scores(pc.data), pch = 16, col = as.numeric(sex),
       main = "Normalized Liver Metabolites Colored by Sex")
  legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
  
  # layout(matrix(1:2, 1, 2))
  
  batch = factor(data.out$Batch)
  plot(scores(pc.data), pch = 16, col = as.numeric(batch),
       main = "Normalized Liver Metabolites Colored by Batch")
  legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
         x.intersp = 0.7, y.intersp = 0.7)
  
  # layout(matrix(1:2, 1, 2))
  wave = factor(data.out$Wave)
  plot(scores(pc.data), pch = 16, col = as.numeric(wave),
       main = "Normalized Liver Metabolites Colored by Wave")
  legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
         x.intersp = 0.7, y.intersp = 0.7)
  
  # layout(matrix(1:2, 1, 2))
  # diet.days = factor(data.out$diet_days, levels = sort(unique(data.out$diet_days)))
  # diet.colors = rainbow(length(levels(diet.days)) - 1)
  # plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
  #     main = "Normalized Liver Metabolites Colored by Diet Days")
  # plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
  #     main = "Normalized Liver Metabolites Colored by Diet Days")

dev.off()

##############################################################################################
# Look at the distribution of phenotypes and the correlation between phenotypes and samples. #
##############################################################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

annot.data.out = data.out[,1:6]
data.data.out  = log2(as.matrix(data.out[,-(1:6)]))

pdf("figures/liver_metabolites_normalized_boxplot.pdf", width = 12, height = 7)

  boxplot(data.data.out, range = 0)
  
dev.off()

pdf("figures/liver_metabolites_normalized_heatmap.pdf", width = 12, height = 12)

  batch.colors = rainbow(12)[as.numeric(factor(annot$Batch))]
  heatmap(data.data.out, RowSideColors = batch.colors)
  
dev.off()


pdf("figures/liver_metabolites_normalized_heatmap.pdf", width = 12, height = 12)

  batch.colors = rainbow(12)[as.numeric(factor(annot$Batch))]
  heatmap(data.data.out, RowSideColors = batch.colors)
  
dev.off()
