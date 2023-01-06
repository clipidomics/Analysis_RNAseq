### PATHWAY ANALYSIS #####
##########################

 library(DESeq2)
 library(clusterProfiler)
 library(biomaRt)
 library(ReactomePA)
 library(DOSE)
 library(KEGG.db)
 library(org.Mm.eg.db)
# library(org.Hs.eg.db)
 library(pheatmap)
 library(genefilter)
 library(RColorBrewer)
 library(GO.db)
 library(topGO)
 library(dplyr)
#BiocManager::install("gage") ; 
 library(gage)
 library(ggsci)

#Mouse genome database (Select the correct one)
results_annotated<-sig_res
rownames(results_annotated)<-results_annotated$gene
# Add gene full name
results_annotated$description <- mapIds(x = org.Mm.eg.db,
                              keys = rownames(results_annotated),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

# Add gene symbol
#results_annotated$symbol <- results_annotated$gene
rownames(results_annotated)<-results_annotated$gene
# Add ENTREZ ID
results_annotated$entrez <- mapIds(x = org.Mm.eg.db,
                         keys = rownames(results_annotated),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

# Add ENSEMBL
rownames(results_annotated)<-results_annotated$gene
results_annotated$ensembl <- mapIds(x = org.Mm.eg.db,
                          keys = row.names(results_annotated),
                          column = "ENSEMBL",
                          keytype = "SYMBOL",
                          multiVals = "first")
#Remove any genes that do not have any entrez identifiers

results_annotated_sig_entrez <- subset(results_annotated, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- results_annotated_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results_annotated_sig_entrez$entrez

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)
library(clusterProfiler)

#KEGG enrichment
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'mouse',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

## Go to see terms kegg_enrich@result

# Plot results_annotated
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)


#GO enrichmentment
go_enrich <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Mm.eg.db', 
                      readable = T,
                      ont = "MF",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

# Plot results_annotated
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 40, 
        title = "GO Biological Pathways",
        font.size = 6)


# Load pathview
library(pathview)

# Plot specific KEGG pathways (with fold change) 
## pathway.id : KEGG pathway identifier
pathview(gene.data = gene_matrix, 
         pathway.id = "00604", 
         species = "mouse")



