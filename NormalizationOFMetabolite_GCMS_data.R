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
metab = right_join(annot, DOLiver01Aug2017, by = "Mouse.ID") # Merge the sample annotation and data by mouse id.
x <- metab[grepl("Redo", metab$Mouse.ID),] # 14 samples with REDO in identifier
metab <- metab[!grepl("Redo", metab$Mouse.ID),] # removing all the rows with Redo in identifier
row.names(metab) <- metab[,1] # set row names
dim(metab) # 384 x 385

# annot = as.data.frame(metab[,1:24])
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
str(pc.data)

pdf("Figures/Liver_metabolites_Tier4_Log2RAW_DOMice_PCA_20170802_2.pdf") # save pdf in Figures folder
  
  batch.colors = as.numeric(factor(metab$Batch.x))
  plot(scores(pc.data), pch = 16, col = batch.colors, main = "Un-normalized Liver Metabolites, Colored by Batch 20170802")
  text(scores(pc.data)[,1], scores(pc.data)[,2], labels = rownames(data.log2), 
       col = batch.colors)

dev.off()

prop.missing = rowMeans(is.na(data)) # mean of each row where there is missing data
sum(prop.missing > 0.25) # determin what samples have more than 25% of missing data
rownames(data)[prop.missing > 0.25] # [1] "DO109" "DO134"
keep = which(prop.missing < 0.25) # keep all data that has less than 25% of missing data
data = data[keep,] # subset data to keep
annot = annot[keep,] # subset annotations to keep
#annot$sex[annot$sex == "XO"] = "F" 

pdf("Figures/VariationbyLoading_LiverMetabolites_20170803.pdf", useDingbats = FALSE) # save pdf in Figures folder
  pc.data = pca(data.log2, method = "bpca", nPcs = 20)
  plot(pc.data)
  abline(h = 0.95, col = 2)
dev.off()

pdf("Figures/LoadingUnnomarlizedLiverMetabolites20170803.pdf", useDingbats = FALSE) # save pdf in Figures folder
  plot(loadings(pc.data), pch = 19, col = "mediumorchid2", main = "Loadings Un-normalized Liver Metabolites, 363 Metabolites") #plot the loadings
  text(loadings(pc.data)[,1], loadings(pc.data)[,2], labels = colnames(data.log2), 
       col = "midnightblue")
dev.off()
dim(loadings(pc.data)) #363 x 20

######################################################
# PCA plots of Unnormalized data by batch, sex, etc. #
######################################################

sex = factor(metab$sex) # create a factor for sex

pdf("figures/liver_metabolites_unnormalized_PCA_20170804.pdf", useDingbats = FALSE) # compiles the following plot into one pdf

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
mod = model.matrix(~sex, data = annot)[,-3]
batch = annot$Batch

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
  data.cb = ComBat(dat = t(data.compl), batch = batch, mod = mod)
  data.cb = t(data.cb)

  # Calculate the change.
  chg = sum((data.compl[miss] - data.cb[miss])^2)
  print(paste("   SS Change:", chg))

  # Put the missing data back in an impute again.
  if(chg > 1 & iter < 20) {

    data.cb[miss] = NA
    data.log = data.cb
    iter = iter + 1    

  } else {

    data.log = data.cb
    break
  }}
)

# Remove or average duplicate samples.
dupl = which(duplicated(rownames(data.log)))
dupl.data = data.log[rownames(data.log) %in% rownames(data.log)[dupl],]

# Keep the sample with the lowest no-call rate.
stopifnot(rownames(data) == rownames(data.log))
prop.missing = rowMeans(is.na(data))
unique.samples = unique(rownames(data.log))
keep = rep(FALSE, nrow(data.log))

for(i in 1:length(unique.samples)) {

  sample = unique.samples[i]
  wh = which(rownames(data.log) == sample)
  wh = wh[which.min(prop.missing[wh])]
  keep[wh] = TRUE

} # for(i)

data.log = data.log[keep,]
annot = annot[match(rownames(data.log), annot$Mouse.ID),]

# Merge in the Chr M and Y info.
attie_MY = read_csv(paste0(input.dir, "attie_sample_info_ChrM_Y.csv"))
annot = right_join(annot, attie_MY, by = "Mouse.ID")
annot = annot[,c(1:10, 13:15)]
colnames(annot) = sub("\\.x", "", colnames(annot))

data.log = data.frame(Mouse.ID = rownames(data.log), data.log)
data.out = right_join(annot, data.log, by = "Mouse.ID")

saveRDS(data.out, file = paste0(output.dir, "attie_liver_metabolites_normalized.rds"))


#################################################
# Do this step if want rank Z-transformed data. # 
#################################################

# Transform each analyte into Z-scores.
data.rz = data.out

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

for(i in 14:ncol(data.rz)) {
  data.rz[,i] = rankZ(data.rz[,i])
}

saveRDS(data.rz, file = paste0(output.dir, "attie_liver_metabolites_zscore_normalized.rds"))

####################################################
# PCA plots of Normalized data by batch, sex, etc. #
####################################################

pdf("figures/liver_metabolites_normalized_PCA.pdf", width = 12, height = 7)

  pc.data = pca(as.matrix(data.out[,-(1:13)]), method = "bpca", nPcs = 20)
  
  layout(matrix(1:2, 1, 2))
  sex = factor(data.out$sex)
  plot(scores(pc.data), pch = 16, col = as.numeric(sex),
       main = "Normalized Liver Metabolites Colored by Sex")
  legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
       main = "Normalized Liver Metabolites Colored by Sex")
  legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
  
  layout(matrix(1:2, 1, 2))
  batch = factor(data.out$Batch)
  plot(scores(pc.data), pch = 16, col = as.numeric(batch),
       main = "Normalized Liver Metabolites Colored by Batch")
  legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
         x.intersp = 0.7, y.intersp = 0.7)
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
       main = "Normalized Liver Metabolites Colored by Batch")
  legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
         x.intersp = 0.7, y.intersp = 0.7)
  
  layout(matrix(1:2, 1, 2))
  wave = factor(data.out$wave)
  plot(scores(pc.data), pch = 16, col = as.numeric(wave),
       main = "Normalized Liver Metabolites Colored by Wave")
  legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
         x.intersp = 0.7, y.intersp = 0.7)
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
       main = "Normalized Liver Metabolites Colored by Wave")
  legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
         x.intersp = 0.7, y.intersp = 0.7)
  
  layout(matrix(1:2, 1, 2))
  diet.days = factor(data.out$diet_days, levels = sort(unique(data.out$diet_days)))
  diet.colors = rainbow(length(levels(diet.days)) - 1)
  plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
       main = "Normalized Liver Metabolites Colored by Diet Days")
  plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
       main = "Normalized Liver Metabolites Colored by Diet Days")

dev.off()

##############################################################################################
# Look at the distribution of phenotypes and the correlation between phenotypes and samples. #
##############################################################################################

annot = data.out[,1:13]
data  = as.matrix(data.out[,-(1:13)])

pdf("figures/liver_metabolites_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(data, range = 0)
dev.off()

pdf("figures/liver_metabolites_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot$Batch))]
heatmap(data, RowSideColors = batch.colors)
dev.off()

################################################################################
# Compare to U. Wisc. data.
norm = read.delim(paste0(input.dir,"14June2017_DOLiverMetabolites_NORM.txt"))
rownames(norm) = norm$Mouse.ID

# Transform each analyte into Z-scores.
norm.rz = norm

for(i in 2:ncol(norm.rz)) {
  norm.rz[,i] = rankZ(norm.rz[,i])
}

saveRDS(norm.rz, file = paste0(output.dir, "attie_liver_metabolites_zscore_uwisc_normalized.rds"))

norm = as.matrix(norm[,-1])

# Merge the sample annotation with the U. Wisc. normalized data.
annot.wisc = annot[annot$Mouse.ID %in% rownames(norm),]
norm = norm[annot.wisc$Mouse.ID,]
stopifnot(annot.wisc$Mouse.ID == rownames(norm))

pc.data = pca(norm, method = "bpca", nPcs = 20)

pdf("figures/liver_metabolites_UWisc_normalized_PCA.pdf", width = 12, height = 7)

layout(matrix(1:2, 1, 2))
sex = factor(annot.wisc$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Liver Metabolites Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Liver Metabolites Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(annot.wisc$Batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Liver Metabolites Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Liver Metabolites Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(annot.wisc$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Liver Metabolites Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Liver Metabolites Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(annot.wisc$diet_days, levels = sort(unique(annot.wisc$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Liver Metabolites Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Liver Metabolites Colored by Diet Days")

dev.off()


# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
pdf("figures/liver_metabolites_UWisc_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(norm, range = 0)
dev.off()

pdf("figures/liver_metabolites_UWisc_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot.wisc$Batch))]
heatmap(norm, RowSideColors = batch.colors)
dev.off()


################################################################################
################################################################################
# Normalize and impute missing data in the Attie liver lipid data set.
# Daniel Gatti
# dan.gatti@jax.org
# July 14, 2017
################################################################################
################################################################################

options(stringsAsFActors = F)
library(tidyverse)
library(pcaMethods)
library(sva)

input.dir = "C:/Users/mkeller3/Desktop/DO_liver/data/"
output.dir = "C:/Users/mkeller3/Desktop/DO_liver/analysis/"
setwd("C:/Users/mkeller3/Desktop/DO_liver/")
getwd()

# Read in the raw lipid data.
lipid = read_delim(paste0(input.dir, "21June2017_DOLiverLipidomicsRawMPK.txt"),
        delim = "\t")

# Read in the sample annotation.
annot = read_delim(paste0(input.dir, "attie_DO_sample_annot.txt"), delim = "\t")

# Merge the sample annotation and data.
lipid = right_join(annot, lipid, by = "Mouse.ID")

# Split up the sample annotation from the data and convert the data into a 
# numeric matrix.
annot = as.data.frame(lipid[,1:11])
data  = as.matrix(lipid[,-(1:11)])
rownames(data)  = annot$Mouse.ID

dim(data)

# One control sample has a 0 value. Set = 1
data[data < 0.01] = 1
range(data)

# Make a PCA plot of all of the data, with sample labels.
pc.data = pca(log(data), method = "bpca", nPcs = 20)

pdf("figures/liver_lipids_unnormalized_all_data_PCA.pdf")

batch.colors = as.numeric(factor(annot$Batch))
plot(scores(pc.data), pch = 16, col = 0, main = "Un-normalized liver Lipids, Colored by Batch")
text(scores(pc.data)[,1], scores(pc.data)[,2], labels = rownames(data), 
     col = batch.colors)

dev.off()

# Remove control samples.
ctrl = which(annot$Mouse.ID == "Control")
data = data[-ctrl,]
annot = annot[-ctrl,]

# 384 samples and 1354 analytes.
dim(data)

######################
# Impute missing data.
data.log = log(data)

# pcaMethods wants samples in rows and variables in columns.
pc.data = pca(data.log, method = "bpca", nPcs = 20)
plot(pc.data)
abline(h = 0.95, col = 2)

# Make PCA plots of the unnormalized data, colored by batch, sex, etc.
pdf("figures/liver_lipids_unnormalized_PCA.pdf")

sex = factor(annot$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Un-normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

batch = factor(annot$Batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "Un-normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       y.intersp = 0.7)

wave = factor(annot$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Un-normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)))

diet.days = factor(annot$diet_days, levels = sort(unique(annot$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Un-normalized Liver Lipids Colored by Diet Days")

dev.off()

boxplot(scores(pc.data)[,1] ~ batch)

# Set up batch and model for comBat.
mod = model.matrix(~sex, data = annot)[,-3]
batch = annot$Batch

# Batch adjust.
# ComBat wants the data with variable in rows and samples in columns.
data.cb = ComBat(dat = t(data.log), batch = batch, mod = mod, prior.plots = TRUE)
data.cb = t(data.cb)

# No duplicate samples.
dupl = which(duplicated(rownames(data.cb)))

# Merge in the Chr M and Y info.
attie_MY = read_csv(paste0(input.dir, "attie_sample_info_ChrM_Y.csv"))
annot = right_join(annot, attie_MY, by = "Mouse.ID")
annot = annot[,c(1:10, 13:15)]
colnames(annot) = sub("\\.x", "", colnames(annot))

data.cb  = data.frame(Mouse.ID = rownames(data.cb), data.cb)
data.out = right_join(annot, data.cb, by = "Mouse.ID")

saveRDS(data.out, file = paste0(output.dir, "attie_liver_lipids_normalized.rds"))

# Transform each analyte into Z-scores.
data.rz = data.out

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

for(i in 14:ncol(data.rz)) {
  data.rz[,i] = rankZ(data.rz[,i])
}

saveRDS(data.rz, file = paste0(output.dir, "attie_liver_lipids_zscore_normalized.rds"))


# Make PCA plots of the normalized data, colored by batch, sex, etc.
pdf("figures/liver_lipids_normalized_PCA.pdf", width = 12, height = 7)

pc.data = pca(as.matrix(data.out[,-(1:13)]), method = "bpca", nPcs = 20)

layout(matrix(1:2, 1, 2))
sex = factor(data.out$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "Normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(data.out$Batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "Normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
     main = "Normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(data.out$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "Normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(data.out$diet_days, levels = sort(unique(data.out$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Normalized Liver Lipids Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "Normalized Liver Lipids Colored by Diet Days")

dev.off()

# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
annot = data.out[,1:13]
data  = as.matrix(data.out[,-(1:13)])

pdf("figures/liver_lipids_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(data, range = 0)
dev.off()

pdf("figures/liver_lipids_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot$Batch))]
heatmap(data, RowSideColors = batch.colors)
dev.off()


###########################
# Compare to U. Wisc. data.
norm = read.delim(paste0(input.dir, "15June2017_DOLiverLipidomics.txt"))
rownames(norm) = norm$Mouse.ID

# Transform each analyte into Z-scores.
norm.rz = norm

for(i in 2:ncol(norm.rz)) {
  norm.rz[,i] = rankZ(norm.rz[,i])
}

saveRDS(norm.rz, file = paste0(output.dir, "attie_liver_lipids_zscore_uwisc_normalized.rds"))

norm = as.matrix(norm[,-1])


# Merge the sample annotation with the U. Wisc. normalized data.
annot.wisc = annot[annot$Mouse.ID %in% rownames(norm),]
norm = norm[annot.wisc$Mouse.ID,]
stopifnot(annot.wisc$Mouse.ID == rownames(norm))

pc.data = pca(norm, method = "bpca", nPcs = 20)

pdf("figures/liver_lipids_UWisc_normalized_PCA.pdf", width = 12, height = 7)

layout(matrix(1:2, 1, 2))
sex = factor(annot.wisc$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Liver Lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(annot.wisc$Batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Liver Lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(annot.wisc$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Liver Lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(annot.wisc$diet_days, levels = sort(unique(annot.wisc$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Liver Lipids Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Liver Lipids Colored by Diet Days")

dev.off()

# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
pdf("figures/liver_lipids_UWisc_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(norm, range = 0)
dev.off()

pdf("figures/Liver_lipids_UWisc_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot.wisc$Batch))]
heatmap(norm, RowSideColors = batch.colors)
dev.off()


################################################################################
################################################################################
# Normalize and impute missing data in the Attie plasma lipid data set.
# Daniel Gatti
# dan.gatti@jax.org
# July 14, 2017
################################################################################
################################################################################

options(stringsAsFActors = F)
library(tidyverse)
library(pcaMethods)
library(sva)

input.dir = "C:/Users/mkeller3/Desktop/DO_liver/data/"
output.dir = "C:/Users/mkeller3/Desktop/DO_liver/analysis/"
setwd("C:/Users/mkeller3/Desktop/DO_liver/")
getwd()

# Read in the raw lipid data.
lipid = read_delim(paste0(input.dir, "21June2017_DOPlasmaLipidomicsRaw.txt"),
        delim = "\t")

# Read in the sample annotation.
annot = read_delim(paste0(input.dir, "attie_DO_sample_annot.txt"), delim = "\t")

# Merge the sample annotation and data.
lipid = right_join(annot, lipid, by = "Mouse.ID")

# Split up the sample annotation from the data and convert the data into a 
# numeric matrix.
annot = as.data.frame(lipid[,1:11])
data  = as.matrix(lipid[,-(1:11)])
rownames(data)  = annot$Mouse.ID

dim(data)

# identify values that are zero
zero_data=which(data == 0)
data[zero_data]=1
zero_data

# Make a PCA plot of all of the data, with sample labels.
pc.data = pca(log(data), method = "bpca", nPcs = 20)

pdf("figures/Plasma_Lipids_Unnormalized_All_Data_PCA.pdf")

batch.colors = as.numeric(factor(annot$Batch))
plot(scores(pc.data), pch = 16, col = 0, main = "Un-normalized Plasma Lipids, Colored by Batch")
text(scores(pc.data)[,1], scores(pc.data)[,2], labels = rownames(data), 
     col = batch.colors)

dev.off()

# Remove control samples.
ctrl = which(annot$Mouse.ID == "Control")
data = data[-ctrl,]
annot = annot[-ctrl,]

# 384 samples and 1354 analytes.
dim(data)

######################
# Impute missing data.
data.log = log(data)

# pcaMethods wants samples in rows and variables in columns.
pc.data = pca(data.log, method = "bpca", nPcs = 20)
plot(pc.data)
abline(h = 0.95, col = 2)

# Make PCA plots of the unnormalized data, colored by batch, sex, etc.
pdf("figures/plasma_lipids_unnormalized_PCA.pdf")

sex = factor(annot$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Un-normalized plasma lipids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

batch = factor(annot$Batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "Un-normalized plasma lipids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       y.intersp = 0.7)

wave = factor(annot$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Un-normalized plasma lipids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)))

diet.days = factor(annot$diet_days, levels = sort(unique(annot$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Un-normalized plasma lipids Colored by Diet Days")

dev.off()

boxplot(scores(pc.data)[,1] ~ batch)

# Set up batch and model for comBat.
mod = model.matrix(~sex, data = annot)[,-3]
batch = annot$Batch

# Batch adjust.
# ComBat wants the data with variable in rows and samples in columns.
data.cb = ComBat(dat = t(data.log), batch = batch, mod = mod, prior.plots = TRUE)
data.cb = t(data.cb)

# No duplicate samples.
dupl = which(duplicated(rownames(data.cb)))

# Merge in the Chr M and Y info.
attie_MY = read_csv(paste0(input.dir, "attie_sample_info_ChrM_Y.csv"))
annot = right_join(annot, attie_MY, by = "Mouse.ID")
annot = annot[,c(1:10, 13:15)]
colnames(annot) = sub("\\.x", "", colnames(annot))

data.cb  = data.frame(Mouse.ID = rownames(data.cb), data.cb)
data.out = right_join(annot, data.cb, by = "Mouse.ID")

saveRDS(data.out, file = paste0(output.dir, "attie_plasma_lipids_normalized.rds"))

# Transform each analyte into Z-scores.
data.rz = data.out

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

for(i in 14:ncol(data.rz)) {
  data.rz[,i] = rankZ(data.rz[,i])
}

saveRDS(data.rz, file = paste0(output.dir, "attie_plasma_lipids_zscore_normalized.rds"))


# Make PCA plots of the normalized data, colored by batch, sex, etc.
pdf("figures/plasma_lipids_normalized_PCA.pdf", width = 12, height = 7)

pc.data = pca(as.matrix(data.out[,-(1:13)]), method = "bpca", nPcs = 20)

  layout(matrix(1:2, 1, 2))
  sex = factor(data.out$sex)
  plot(scores(pc.data), pch = 16, col = as.numeric(sex),
       main = "Normalized Plasma Lipids Colored by Sex")
  legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
       main = "Normalized Plasma Lipids Colored by Sex")
  legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
  
  layout(matrix(1:2, 1, 2))
  batch = factor(data.out$Batch)
  plot(scores(pc.data), pch = 16, col = as.numeric(batch),
       main = "Normalized Plasma Lipids Colored by Batch")
  legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
         x.intersp = 0.7, y.intersp = 0.7)
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
       main = "Normalized Plasma Lipids Colored by Batch")
  legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
         x.intersp = 0.7, y.intersp = 0.7)
  
  layout(matrix(1:2, 1, 2))
  wave = factor(data.out$wave)
  plot(scores(pc.data), pch = 16, col = as.numeric(wave),
       main = "Normalized Plasma Lipids Colored by Wave")
  legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
         x.intersp = 0.7, y.intersp = 0.7)
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
       main = "Normalized Plasma Lipids Colored by Wave")
  legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
         x.intersp = 0.7, y.intersp = 0.7)
  
  layout(matrix(1:2, 1, 2))
  diet.days = factor(data.out$diet_days, levels = sort(unique(data.out$diet_days)))
  diet.colors = rainbow(length(levels(diet.days)) - 1)
  plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
       main = "Normalized Plasma Lipids Colored by Diet Days")
  plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
       main = "Normalized Plasma Lipids Colored by Diet Days")

dev.off()

# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
annot = data.out[,1:13]
data  = as.matrix(data.out[,-(1:13)])

pdf("figures/plasma_lipids_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(data, range = 0)
dev.off()

pdf("figures/plasma_lipids_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot$Batch))]
heatmap(data, RowSideColors = batch.colors)
dev.off()


###########################
# Compare to U. Wisc. data.
norm = read.delim(paste0(input.dir, "15June2017_DOPlasmaLipidomics.txt"))
rownames(norm) = norm$Mouse.ID

# Transform each analyte into Z-scores.
norm.rz = norm

for(i in 2:ncol(norm.rz)) {
  norm.rz[,i] = rankZ(norm.rz[,i])
}

saveRDS(norm.rz, file = paste0(output.dir, "attie_plasma_lipids_zscore_uwisc_normalized.rds"))

norm = as.matrix(norm[,-1])

# Merge the sample annotation with the U. Wisc. normalized data.
annot.wisc = annot[annot$Mouse.ID %in% rownames(norm),]
norm = norm[annot.wisc$Mouse.ID,]
stopifnot(annot.wisc$Mouse.ID == rownames(norm))

pc.data = pca(norm, method = "bpca", nPcs = 20)

pdf("figures/plasma_lipids_UWisc_normalized_PCA_2.pdf", width = 12, height = 7)

  layout(matrix(1:2, 1, 2))
  sex = factor(annot.wisc$sex)
  plot(scores(pc.data), pch = 16, col = as.numeric(sex),
       main = "U. Wisc. Normalized Plasma Lipids Colored by Sex")
  legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
       main = "U. Wisc. Normalized Plasma Lipids Colored by Sex")
  legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
  
  layout(matrix(1:2, 1, 2))
  batch = factor(annot.wisc$Batch)
  plot(scores(pc.data), pch = 16, col = as.numeric(batch),
       main = "U. Wisc. Normalized Plasma Lipids Colored by Batch")
  legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
         x.intersp = 0.7, y.intersp = 0.7)
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
       main = "U. Wisc. Normalized Plasma Lipids Colored by Batch")
  legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
         x.intersp = 0.7, y.intersp = 0.7)
  
  layout(matrix(1:2, 1, 2))
  wave = factor(annot.wisc$wave)
  plot(scores(pc.data), pch = 16, col = as.numeric(wave),
       main = "U. Wisc. Normalized Plasma Lipids Colored by Wave")
  legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
         x.intersp = 0.7, y.intersp = 0.7)
  plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
       main = "U. Wisc. Normalized Plasma Lipids Colored by Wave")
  legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
         x.intersp = 0.7, y.intersp = 0.7)
  
  layout(matrix(1:2, 1, 2))
  diet.days = factor(annot.wisc$diet_days, levels = sort(unique(annot.wisc$diet_days)))
  diet.colors = rainbow(length(levels(diet.days)) - 1)
  plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
       main = "U. Wisc. Normalized Plasma Lipids Colored by Diet Days")
  plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
       main = "U. Wisc. Normalized Plasma Lipids Colored by Diet Days")

dev.off()

# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
pdf("figures/plasma_lipids_UWisc_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(norm, range = 0)
dev.off()

pdf("figures/Plasma_lipids_UWisc_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot.wisc$Batch))]
heatmap(norm, RowSideColors = batch.colors)
dev.off()



################################################################################
################################################################################
# Normalize and impute missing data in the Attie cecum bile acid data set.
# Daniel Gatti
# dan.gatti@jax.org
# July 14, 2017
################################################################################
################################################################################

options(stringsAsFActors = F)
library(tidyverse)
library(pcaMethods)
library(sva)

input.dir = "C:/Users/mkeller3/Desktop/DO_liver/data/"
output.dir = "C:/Users/mkeller3/Desktop/DO_liver/analysis/"
setwd("C:/Users/mkeller3/Desktop/DO_liver/")
getwd()

# Read in the raw lipid data.
BA = read_delim(paste0(input.dir, "16_June_2017_DO_Cecum_Bile_Acids_RAW.txt"),
        delim = "\t")

# Read in the sample annotation.
annot = read_delim(paste0(input.dir, "attie_DO_sample_annot.txt"), delim = "\t")

# Merge the sample annotation and data.
BA = right_join(annot, BA, by = "Mouse.ID")

colnames(BA)

# Split up the sample annotation from the data and convert the data into a 
# numeric matrix.
annot = as.data.frame(BA[,1:9])
data  = as.matrix(BA[,-(1:9)])
rownames(data)  = annot$Mouse.ID

dim(data)
data[1:5,1:5]
colnames(data)


# identify values that are zero
zero_data=which(data == 0)
data[zero_data]=1
zero_data

# Make a PCA plot of all of the data, with sample labels.
pc.data = pca(log(data), method = "bpca", nPcs = 20)

pdf("figures/Cecum_Bile_Acids_Unnormalized_All_Data_PCA.pdf")

batch.colors = as.numeric(factor(annot$Batch))
plot(scores(pc.data), pch = 16, col = 0, main = "Un-normalized Cecum Bile Acids, Colored by Batch")
text(scores(pc.data)[,1], scores(pc.data)[,2], labels = rownames(data), 
     col = batch.colors)

dev.off()

colnames(data)
plot(data[,29])

#Remove spike-in control samples
spike.in = grep("^CONT ", colnames(data))
colnames(data)
data=data[ ,-spike.in]


# 383 samples and 27 analytes.
dim(data)

######################
# Impute missing data.
data.log = log(data+1)
data.log[1:5,1:5]
hist(data.log)


# pcaMethods wants samples in rows and variables in columns.
pc.data = pca(data.log, method = "bpca", nPcs = 20)
plot(pc.data)
abline(h = 0.95, col = 2)

# Make PCA plots of the unnormalized data, colored by batch, sex, etc.
pdf("figures/Cecum_Bile_Acids_unnormalized_PCA.pdf")

sex = factor(annot$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Un-normalized Cecum_Bile_Acids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

batch = factor(annot$Batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "Un-normalized Cecum_Bile_Acids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       y.intersp = 0.7)

wave = factor(annot$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Un-normalized Cecum_Bile_Acids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)))

diet.days = factor(annot$diet_days, levels = sort(unique(annot$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Un-normalized Cecum_Bile_Acids Colored by Diet Days")

dev.off()

boxplot(scores(pc.data)[,1] ~ batch)

# Set up batch and model for comBat.
colnames(annot)


# Set up batch and model for comBat.
mod = model.matrix(~sex, data = annot)[,-3]
batch = annot$Batch

# Imputing missing values and batch normalization
chg = 1e6
iter = 1
repeat( {

  print(paste("Iteration", iter))

  # Impute missing data.
  miss = which(is.na(data.log))
  print(paste(length(miss), "missing points."))
  pc.data = pca(data.log, method = "bpca", nPcs = 7)
  data.compl = completeObs(pc.data)

  # Batch adjust.
  # ComBat wants the data with variable in rows and samples in columns.
  data.cb = ComBat(dat = t(data.compl), batch = batch, mod = mod)
  data.cb = t(data.cb)

  # Calculate the change.
  chg = sum((data.compl[miss] - data.cb[miss])^2)
  print(paste("   SS Change:", chg))

  # Put the missing data back in an impute again.
  if(chg > 1 & iter < 20) {

    data.cb[miss] = NA
    data.log = data.cb
    iter = iter + 1    

  } else {

    data.log = data.cb
    break
  }}
)

colnames(data.log)

# No duplicate samples.
dupl = which(duplicated(rownames(data.cb)))
dupl

# Merge in the Chr M and Y info.
attie_MY = read_csv(paste0(input.dir, "attie_sample_info_ChrM_Y.csv"))
annot = right_join(annot, attie_MY, by = "Mouse.ID")
annot = annot[,-grep(".y$",colnames(annot))]
colnames(annot)

colnames(annot) = sub("\\.x", "", colnames(annot))

data.cb  = data.frame(Mouse.ID = rownames(data.cb), data.cb)
data.out = right_join(annot, data.cb, by = "Mouse.ID")

saveRDS(data.out, file = paste0(output.dir, "attie_cecum_bile_acids_normalized.rds"))

# Transform each analyte into Z-scores.
data.rz = data.out

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

for(i in 14:ncol(data.rz)) {
  data.rz[,i] = rankZ(data.rz[,i])
}

saveRDS(data.rz, file = paste0(output.dir, "attie_cecum_bile_acids_zscore_normalized.rds"))


# Make PCA plots of the normalized data, colored by batch, sex, etc.
pdf("figures/cecum_bile_acids_normalized_PCA.pdf", width = 12, height = 7)

pc.data = pca(as.matrix(data.out[,-(1:13)]), method = "bpca", nPcs = 20)

layout(matrix(1:2, 1, 2))
sex = factor(data.out$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "Normalized Cecum Bile Acids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "Normalized Cecum Bile Acids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(data.out$Batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "Normalized Cecum Bile Acids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
     main = "Normalized Cecum Bile Acids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(data.out$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "Normalized Cecum Bile Acids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "Normalized Cecum Bile Acids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(data.out$diet_days, levels = sort(unique(data.out$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "Normalized Cecum Bile Acids Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "Normalized Cecum Bile Acids Colored by Diet Days")

dev.off()

# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
annot = data.out[,1:13]
data  = as.matrix(data.out[,-(1:13)])

pdf("figures/Cecum_Bile_Acids_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(data, range = 0)
dev.off()

pdf("figures/Cecum_Bile_Acids_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot$Batch))]
heatmap(data, RowSideColors = batch.colors)
dev.off()


###########################
# Compare to U. Wisc. data.
norm = read.delim(paste0(input.dir, "16_June_2017_DO_Cecum_Bile_Acids_UW_Norm.txt"))
rownames(norm) = norm$Mouse.ID
colnames(norm)

# Transform each analyte into Z-scores.
norm.rz = norm

for(i in 2:ncol(norm.rz)) {
  norm.rz[,i] = rankZ(norm.rz[,i])
}

saveRDS(norm.rz, file = paste0(output.dir, "attie_cecum_bile_acids_zscore_uwisc_normalized.rds"))

norm = as.matrix(norm[,-(1:3)])
colnames(norm)

# Merge the sample annotation with the U. Wisc. normalized data.
annot.wisc = annot[annot$Mouse.ID %in% rownames(norm),]
norm = norm[annot.wisc$Mouse.ID,]
stopifnot(annot.wisc$Mouse.ID == rownames(norm))

pc.data = pca(norm, method = "bpca", nPcs = 20)

pdf("figures/cecum_bile_acids_UWisc_normalized_PCA.pdf", width = 12, height = 7)

layout(matrix(1:2, 1, 2))
sex = factor(annot.wisc$sex)
plot(scores(pc.data), pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Cecum Bile Acids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(sex),
     main = "U. Wisc. Normalized Cecum Bile Acids Colored by Sex")
legend("bottomleft", legend = levels(sex), pch = 16, col = 1:length(levels(sex)))

layout(matrix(1:2, 1, 2))
batch = factor(annot.wisc$Batch)
plot(scores(pc.data), pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Cecum Bile Acids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(batch),
     main = "U. Wisc. Normalized Cecum Bile Acids Colored by Batch")
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
wave = factor(annot.wisc$wave)
plot(scores(pc.data), pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Cecum Bile Acids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)
plot(scores(pc.data)[,3:2], pch = 16, col = as.numeric(wave),
     main = "U. Wisc. Normalized Cecum Bile Acids Colored by Wave")
legend("bottomleft", legend = levels(wave), pch = 16, col = 1:length(levels(wave)),
       x.intersp = 0.7, y.intersp = 0.7)

layout(matrix(1:2, 1, 2))
diet.days = factor(annot.wisc$diet_days, levels = sort(unique(annot.wisc$diet_days)))
diet.colors = rainbow(length(levels(diet.days)) - 1)
plot(scores(pc.data), pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Cecum Bile Acids Colored by Diet Days")
plot(scores(pc.data)[,3:2], pch = 16, col = diet.colors[diet.days],
     main = "U. Wisc. Normalized Cecum Bile Acids Colored by Diet Days")

dev.off()

# Look at the distribution of phenotypes and the correlation between phenotypes
# and samples.
pdf("figures/cecum_bile_acids_UWisc_normalized_boxplot.pdf", width = 12, height = 7)
boxplot(norm, range = 0)
dev.off()

pdf("figures/Cecum_Bile_acids_UWisc_normalized_heatmap.pdf", width = 12, height = 12)
batch.colors = rainbow(12)[as.numeric(factor(annot.wisc$Batch))]
heatmap(norm, RowSideColors = batch.colors)
dev.off()


####################################################################################
# stop here
# next step will be to add script for plasma BA measurements
####################################################################################










