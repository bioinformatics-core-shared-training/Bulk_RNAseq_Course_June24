library(DESeq2)
library(tidyverse)

txi <- readRDS("RObjects/txi.rds")
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv")

# Check that the column names in the txi match the 
# SampleNames in the samplinfo

all(colnames(txi$counts)==sampleinfo$SampleName)
# TRUE

# Create the model design as a formula

simple.model <-  as.formula(~ Status)

model.matrix(simple.model, data = sampleinfo)

sampleinfo <- mutate(sampleinfo,
    Status = fct_relevel(Status, "Uninfected"))

model.matrix(simple.model, data = sampleinfo)

# Build a DESeq2DataSet object

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)

# Filter out unexpressed (uninformative) genes

counts(ddsObj.raw)

keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep, ]

# Differential Expression Analysis

# Estimate size factors

ddsObj <- estimateSizeFactors(ddsObj.filt)

## Looking at the effect of size factors

logcounts <- log2(counts(ddsObj, normalized = FALSE) + 1)
logNormCounts <-log2(counts(ddsObj, normalized = TRUE) + 1)

limma::plotMA(logcounts, array = 5, ylim = c(-5, 5))
abline(h = 0, col = "red")


limma::plotMA(logNormCounts, array = 5, ylim = c(-5, 5))
abline(h = 0, col = "red")

# Estimate Dispersion 

ddsObj <- estimateDispersions(ddsObj)
plotDispEsts(ddsObj)

# Fit model and apply the statistical test 

## Wald Test

ddsObj <- nbinomWaldTest(ddsObj)

# The `DESeq` command

ddsObj <- DESeq(ddsObj.filt)

# Get the results of the DE analysis

results.simple <- results(ddsObj, alpha = 0.05)

# Exercise 1

# How many up regulated genes

sum(results.simple$log2FoldChange > 0 & 
  results.simple$padj < 0.05, 
  na.rm = TRUE)

# How many down regulated genes

sum(results.simple$log2FoldChange < 0 & 
      results.simple$padj < 0.05, 
    na.rm = TRUE)

# Independent Filtering

sum(is.na(results.simple$padj))
sum(is.na(results.simple$pvalue))

# Additive Model

## Exercise 2

additive.model <- as.formula( ~ TimePoint + Status )

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
        colData = sampleinfo,
        design = additive.model)

ddsObj.filt <- ddsObj.raw[keep, ]
ddsObj <- DESeq(ddsObj.filt)

# Extract the restults 

results.additive <- results(ddsObj, alpha = 0.05)

sum(results.additive$padj < 0.05, na.rm =TRUE)

hist(results.additive$pvalue)

# Retrieving results for other contrasts

results.InfectedvUninfected <- results.additive
rm(results.additive)

resultsNames(ddsObj)

# Exercise 3

# retrieve results for TimePoint d33 v d11

resultsNames(ddsObj)

results.d33vd11 <- results(ddsObj,
     alpha = 0.05,
     name = "TimePoint_d33_vs_d11")


sum(results.d33vd11$padj < 0.05, na.rm = TRUE)

## Comparing two design models

ddsObj.LRT <- DESeq(ddsObj, 
                    test = "LRT", 
                    reduced = simple.model)

results.Additive_v_Simple <- 
  results(ddsObj.LRT, alpha = 0.05)

# Exercise 4

## Create new filtered DESeqDataSet Object

interaction.model <- 
  as.formula(~ TimePoint + Status + TimePoint:Status)

### Modify the design formula
design(ddsObj.filt) <- interaction.model

### Run DESeq2
ddsObj.interaction <- DESeq(ddsObj.filt)

## 3. Use the LRT

ddsObj.LRT <- DESeq(ddsObj.interaction, 
                    test = "LRT",
                    reduced = additive.model)

results.Interaction_v_Additive <- results(ddsObj.LRT,
                                          alpha = 0.05)

### 4. How many significant

table(results.Interaction_v_Additive$padj < 0.05)

# Extracting the results from the Interaction Model

resultsNames(ddsObj.interaction)

results.interaction.11 <- results(ddsObj.interaction,
         name = "Status_Infected_vs_Uninfected",
         alpha  = 0.05)

model.matrix(interaction.model, data = sampleinfo)

results.interaction.33 <- results(ddsObj.interaction,
  contrast = list(c("Status_Infected_vs_Uninfected",
                    "TimePointd33.StatusInfected")),
  alpha = 0.05)

sum(results.interaction.11$padj < 0.05, na.rm =TRUE)
sum(results.interaction.33$padj < 0.05, na.rm =TRUE)

# Exercise 5

## 1. d33 v d11 for Uninfected mice

results.int.uninf <- results(ddsObj.interaction,
            name = "TimePoint_d33_vs_d11",
            alpha = 0.05)

sum(results.int.uninf$padj < 0.05, na.rm = TRUE)

## 2. d33 v d11 for Infected mice

results.int.inf <- results(ddsObj.interaction,
 contrast = list(c("TimePoint_d33_vs_d11",
                   "TimePointd33.StatusInfected")),
 alpha = 0.05)

sum(results.int.inf$padj < 0.05, na.rm = TRUE)












