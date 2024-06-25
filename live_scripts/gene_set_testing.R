# Gene set testing - ORA and GSEA

# Load libraries ----
library(tidyverse)
library(clusterProfiler)
#for plotting KEGG pathways
library(pathview)
# for mouse GO terms 
library(org.Mm.eg.db)


# 1 - ORA using the KEGG database ----

# Find organism code in KEGG ----
search_kegg_organism('mouse', by = 'common_name')


# Load our DEGs ---- 

shrink.d11 <- readRDS('RObjects/Shrunk_Results.d11.rds')
shrink.d11


# make a vector of Entrex IDs of the DEGs
# p.adj < 0.05 and LFC > 1


sigGenes <- shrink.d11 %>%
  drop_na(Entrez, padj, log2FoldChange) %>%
  filter(padj < 0.05 & abs(log2FoldChange) >1 ) %>%
  pull(Entrez)
  


# Do the ORA 
keggRes <-  enrichKEGG( gene = sigGenes,
                        organism = 'mmu')

keggRes

view(as_tibble(keggRes))

  

# get the pathway in the browser ----
browseKEGG(keggRes, 'mmu04612')


# make a named vector of log fold changes with Entrez IDs
# as the name 
logFC <- shrink.d11$log2FoldChange
names(logFC) <- shrink.d11$Entrez
pathview(gene.data = logFC,
         pathway.id = 'mmu04612',
         species = 'mmu',
         limit = list(gene=20, cpd = 1))

# Exercise 1
# for mmu04659 plot just the sig genes 

logFC_sig <- shrink.d11 %>%
  drop_na(Entrez, padj, log2FoldChange) %>%
  filter(padj < 0.01) %>%
  pull(log2FoldChange, Entrez)
logFC_sig
pathview(gene.data = logFC_sig,
         pathway.id = 'mmu04659',
         species = 'mmu',
         limit = list(gene=20, cpd = 1))


# 2 - ORA using GO (gene ontology) terms

# we need 1 - the gene list of interest 2 - the background (universe)

# 1 gene list
sigGenes_GO <- shrink.d11 %>%
  drop_na(padj, GeneID) %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  pull(GeneID)
sigGenes_GO


# 2 - the universe 
universe <- shrink.d11$GeneID


ego <- enrichGO(gene = sigGenes_GO,
                universe = universe,
                OrgDb = org.Mm.eg.db,
                keyType = 'ENSEMBL',
                ont = 'BP',
                pvalueCutoff = 0.01,
                readable = TRUE)
ego


# make a barplot of the enriched terms
barplot(ego, showCategory = 20)


dotplot(ego, font.size =12)


library(enrichplot)
ego_pt <- pairwise_termsim(ego)
emapplot(ego_pt, cex.params = list(category_label = 0.8))


# 3 - GSEA form the Borad Institute 

library(msigdb) # This contains some code that grabs the gene sets
library(ExperimentHub) # This allow us to pull and do the enricment

# set up the collection
eh <- ExperimentHub()
query(eh, c('msigdb',
            'mm',
            '2023'))


msigdb.mm <- getMsigdb( org = 'mm',
                        id  = 'EZID',
                        version = '2023.1')

msigdb.mm


listCollections(msigdb.mm)


# 1 rank the genes according to some metric
ranked <- shrink.d11 %>%
  drop_na(Entrez, padj, log2FoldChange) %>%
  mutate(rank = log2FoldChange) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)
ranked  
  
  
# 2 load the gene-sets
hallmarks  <- subsetCollection(msigdb.mm, 'h')
hallmarks
msigbd_ids <- geneIds(hallmarks)
msigbd_ids


term2gene <- enframe(msigbd_ids, name = 'gs_name',
                     value = 'entrez') %>%
  unnest(entrez)

term2gene

# 3 do the GSEA

gseaRes <- GSEA(geneList = ranked,
                TERM2GENE = term2gene,
                pvalueCutoff = 1.00,
                minGSSize = 15,
                maxGSSize = 500)

gseaRes


as_tibble(gseaRes) %>%
  arrange(desc(NES)) %>%
  top_n(10, wt=-p.adjust) %>%
  dplyr::select(-core_enrichment)


gseaplot(gseaRes,
         geneSetID = "HALLMARK_INFLAMMATORY_RESPONSE", 
         title = "HALLMARK_INFLAMMATORY_RESPONSE")



# gsea 1 - ranked list - 2 - gene sets 3 - do the analysis

ranked <- shrink.d11 %>%
  drop_na(Entrez, padj, log2FoldChange) %>%
  mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)
ranked  


gseaRes <- GSEA(geneList = ranked,
                TERM2GENE = term2gene,
                pvalueCutoff = 1.00,
                minGSSize = 15,
                maxGSSize = 500)


View(as_tibble(gseaRes) %>%
  arrange(desc(NES)) %>%
  top_n(10, wt=-p.adjust) %>%
  dplyr::select(-core_enrichment))


shrink.d33 <- readRDS('./RObjects/Shrunk_Results.d33.rds')



ranked <- shrink.d33 %>%
  drop_na(Entrez, padj, log2FoldChange) %>%
  mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)
ranked  


gseaRes <- GSEA(geneList = ranked,
                TERM2GENE = term2gene,
                pvalueCutoff = 1.00,
                minGSSize = 15,
                maxGSSize = 500)


View(as_tibble(gseaRes) %>%
       arrange(desc(NES)) %>%
       top_n(10, wt=-p.adjust) %>%
       dplyr::select(-core_enrichment))








