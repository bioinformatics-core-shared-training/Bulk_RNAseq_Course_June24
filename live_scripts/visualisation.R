library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(DESeq2)
library(tidyverse)

# ---- load data ----
ddsObj.interaction <- readRDS("RObjects/DESeqDataSet.interaction.rds")
results.interaction.11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")
results.interaction.33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")

# ---- add annotation to DESeq2 results ----
# create an annotation hub instance
ah <- AnnotationHub()
ah[1]
# download the database we want to use
MouseEnsDb <- query(ah, c("EnsDb", "Mus musculus", "102"))[[1]]

# turn the whole gene-level annotation table into a data frame
annotations <- genes(MouseEnsDb, return.type = "data.frame")
colnames(annotations)

annot <- annotations %>% 
  select(gene_id, gene_name, entrezid) %>% 
  filter(gene_id %in% rownames(results.interaction.11))

sum(is.na(annot$entrezid))

# load pre-curated annotation table
ensemblAnnot <- readRDS("RObjects/Ensembl_annotations.rds")
colnames(ensemblAnnot)

annot.interaction.11 <- as.data.frame(results.interaction.11) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID")

# output the annotated DE results
write_tsv(annot.interaction.11, "results/Interaction.11_Results_Annotated.txt")

# ---- Visualization ----

# p-value histogram
hist(annot.interaction.11$pvalue)

# ---- shrinking log2 fold change ----
ddsShrink.11 <- lfcShrink(ddsObj.interaction,
          res = results.interaction.11,
          type = "ashr")

shrinkTab.11 <- as.data.frame(ddsShrink.11) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID")

# comparison before & after LFC shrinkage
head(results.interaction.11)
head(ddsShrink.11)

# ---- MA plots ----
# show the log fold change for each gene against its average expression
par(mfrow=c(1,2)) # layout will be 1 row, 2 cols
plotMA(results.interaction.11, alpha = 0.05)
plotMA(ddsShrink.11, alpha = 0.05)

# ---- volcano plots ----
ggplot(shrinkTab.11, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), size = 1) +
  geom_text(data = ~top_n(.x, 1, wt = -padj),
            aes(label = Symbol)) +
  labs(x = "log2(fold change)", y = "-log10(p-value)",
       color = "FDR < 5%", title = "Infected vs Uninfected (day 11)")

# ---- Exercise 1 - Volcano plot for 33 days ----
# shrink the results for the 33 days contrast
ddsShrink.33 <- lfcShrink(ddsObj.interaction,
                          res = results.interaction.33,
                          type = "ashr")
shrinkTab.33 <- as.data.frame(ddsShrink.33) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID")

# create a volcano plot
ggplot(shrinkTab.33, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), size = 1) +
  geom_text(data = ~top_n(.x, 5, wt = -padj), aes(label = Symbol)) +
  labs(x = "log2(fold change)", y = "-log10(p-value)",
       color = "FDR < 5%", title = "Infected vs Uninfected (day 33)")

# if interested in down-regulated genes
ggplot(shrinkTab.33, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), size = 1) +
  geom_text(data = shrinkTab.33 %>% filter(log2FoldChange < 0) %>% 
              top_n(2, wt = -padj), aes(label = Symbol)) +
  geom_text(data = ~top_n(.x, 5, wt = -padj), aes(label = Symbol)) +
  labs(x = "log2(fold change)", y = "-log10(p-value)",
       color = "FDR < 5%", title = "Infected vs Uninfected (day 33)")

# ---- MA plot with ggplot2 ----
no_shrink_33 <- ggplot(results.interaction.33,
                       aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), size = .5) +
  scale_color_manual(values = c("TRUE" = "blue",
                                "FALSE" = "grey",
                                "NA" = "grey")) +
  ylim(-4, 4) +
  theme(legend.position = "none")

shrunk_33 <- ggplot(shrinkTab.33,
                    aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), size = .5) +
  scale_color_manual(values = c("TRUE" = "blue",
                                "FALSE" = "grey",
                                "NA" = "grey")) +
  ylim(-4, 4) +
  theme(legend.position = "none")

library(gridExtra)
grid.arrange(no_shrink_33, shrunk_33, ncol = 2)

# ---- strip charts for gene expression ----
# expression levels of individual samples for genes of interest
geneID <- filter(shrinkTab.11, Symbol=="Il10ra") %>% 
  pull(GeneID)

plotCounts(ddsObj.interaction,
           gene = geneID,
           intgroup = c("TimePoint", "Status", "Replicate"),
           returnData = TRUE) %>% 
  ggplot(aes(x=Status, y=log2(count))) +
  geom_point(aes(fill=Replicate), shape=21, size=2) +
  facet_wrap(~TimePoint) +
  expand_limits(y=0) +
  labs(title = "Normalized counts - Interleukin 10 receptor, alpha")

# ---- Exercise 3 - strip chart for the gene 'Jchain' ----
geneID <- filter(shrinkTab.11, Symbol=="Jchain") %>% 
  pull(GeneID)

plotCounts(ddsObj.interaction,
           gene = geneID,
           intgroup = c("TimePoint", "Status", "Replicate"),
           returnData = TRUE) %>% 
  ggplot(aes(x=Status, y=log2(count))) +
  geom_point(aes(fill=Replicate), shape=21, size=2) +
  facet_wrap(~TimePoint) +
  expand_limits(y=0) +
  labs(title = "Normalized counts - Immunoglobulin joining chain")

# ---- Venn Diagram ----
library(ggvenn)
vennDat <- tibble(Geneid = rownames(results.interaction.11)) %>% 
  mutate(Upregulated_11 = results.interaction.11$padj < 0.05 &
           !is.na(results.interaction.11$padj) &
           results.interaction.11$log2FoldChange > 0) %>% 
  mutate(Downregulated_11 = results.interaction.11$padj < 0.05 &
           !is.na(results.interaction.11$padj) &
           results.interaction.11$log2FoldChange < 0) %>% 
  mutate(Upregulated_33 = results.interaction.33$padj < 0.05 &
           !is.na(results.interaction.33$padj) &
           results.interaction.33$log2FoldChange > 0) %>% 
  mutate(Downregulated_33 = results.interaction.33$padj < 0.05 &
           !is.na(results.interaction.33$padj) &
           results.interaction.33$log2FoldChange < 0)

ggvenn(vennDat, set_name_size = 4)

# ---- Heatmap ----
library(ComplexHeatmap)
library(circlize)

# get the top genes
sigGenes <- shrinkTab.11 %>% 
  top_n(300, wt = -padj) %>% 
  pull("GeneID")

plotDat <- vst(ddsObj.interaction)[sigGenes, ] %>% 
  assay()

z.mat <- t(scale(t(plotDat), center = TRUE, scale = TRUE))

# color palette
myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2, 0, 2), myPalette) # color interpolation

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE)

# split the heat map into clusters and add annotation

ha1 = HeatmapAnnotation(
  df = colData(ddsObj.interaction)[ , c("Status",
                                        "TimePoint")],
  col = list(
    Status = c("Uninfected" = "darkgreen",
               "Infected" = "palegreen"),
    TimePoint = c("d11" = "lightblue",
                  "d33" = "darkblue")
  )
)

Heatmap(z.mat, name = "z-score",
        color = myRamp,
        show_row_names = FALSE,
        split = 3, # split the dendrogram by 3
        top_annotation = ha1)

# save data
saveRDS(annot.interaction.11, file = "results/Annotated_Results.d11.rds")
saveRDS(shrinkTab.11, file = "results/Shrunk_Results.d11.rds")
saveRDS(shrinkTab.33, file = "results/Shrunk_Results.d33.rds")

