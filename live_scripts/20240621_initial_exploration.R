# This is for initial exploration
library(tximport)
library(DESeq2)
library(tidyverse)
##############################
## read in meta data of samples
sampleinfo <- read_tsv("data/samplesheet.tsv", col_types = c("cccc"))
arrange(sampleinfo, Status, TimePoint, Replicate)

## read in the count data <- salmon
files <- file.path("salmon", sampleinfo$SampleName, "quant.sf")
files <- set_names(files, sampleinfo$SampleName)
tx2gene <- read_tsv("references/tx2gene.tsv")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

saveRDS(txi, file = "salmon_outputs/txi.rds")

##################### Exercise 1
tpm <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM")

##############################
## Line 7 - 15 can be replaced by L20 here
txi <- readRDS("salmon_outputs/txi.rds")

## prepare count matrix
rawCounts <- round(txi$counts, 0)
dim(rawCounts)

## filter out genes not expressed in more than 5 samples
keep <- rowSums(rawCounts) > 5
table(keep, useNA="always")

filtCounts <- rawCounts[keep, ]
dim(filtCounts)

boxplot(filtCounts, main = "Raw counts", las = 2)

plot(rowMeans(filtCounts), rowSds(filtCounts),
     main = "Raw counts : sd vs mean",
     xlim = c(0, 10000),
     ylim = c(0, 5000))

## data transformation
logcounts <- log2(filtCounts + 1)
boxplot(logcounts, main = "Log2(counts)", las = 2)

# make a color vector
statusCols <- case_when(sampleinfo$Status == "Infected" ~ "navy",
                        sampleinfo$Status == "Uninfected" ~ "gold")

boxplot(logcounts, main = "Log2(Counts)", xlab = "",
        ylab = "Log2(Counts)", las = 2, col = statusCols)

abline(h = median(logcounts), col = "salmon")

plot(rowMeans(logcounts), rowSds(logcounts), main = "Log2 Counts: sd vs mean")

## vst 
vst_counts <- vst(filtCounts)
boxplot(vst_counts, main = "VST (Counts)", xlab = "",
        ylab = "VST(Counts)", las = 2, col = statusCols)

abline(h = median(vst_counts), col = "salmon")
plot(rowMeans(vst_counts), rowSds(vst_counts), main = "Log2 Counts: sd vs mean")

### PCA
library(ggfortify)

pcDat <- prcomp(t(vst_counts))
autoplot(pcDat)

autoplot(pcDat, data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

library(ggrepel)

autoplot(pcDat, data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
  geom_text_repel(aes(x = PC1, y = PC2, label = SampleName),
                  box.padding = 0.8)


sampleinfo <- mutate(sampleinfo, Status=case_when(
  SampleName=="SRR7657882" ~ "Uninfected",
  SampleName=="SRR7657873" ~ "Infected", 
  TRUE ~ Status))

write_tsv(sampleinfo, 
          "results/SampleInfo_Corrected.txt")

autoplot(pcDat,
         data = sampleinfo, 
         colour="Status", 
         shape="TimePoint",
         size=5)

### clustering
library(ggdendro)
hclDat <-  t(vst_counts) %>%
  dist(method = "euclidean") %>%
  hclust()

ggdendrogram(hclDat, rotate=TRUE)

hclDat2 <- hclDat

hclDat2$labels <- str_c(sampleinfo$Status, ":", sampleinfo$TimePoint)
ggdendrogram(hclDat2, rotate=TRUE)
