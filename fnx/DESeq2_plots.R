### Common DESeq Function ###
# Get object DESeq and make some graphics
### 

library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
#biocLite("KEGG.db") ; library(KEGG.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
library(reshape2)
library(data.table)

packageVersion("org.Mm.eg.db")

options(scipen=0.999) ## handle scientific notation


## load dataset

dds<-ddsMat

#Remove some conditions TODO: make options

unique(dds$sample_id)
unique(dds$GroupVar)

#dds<-dds[,which(!dds$sample_id %in% c("ChowAnimal4","ChowAnimal5","ChowAnimal6"))]

## Remove condition which is outlier clearly 
dds<-dds[,which(!dds$sample_id %in% c("DEGS1_5_43_189_004_S41"))]

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
colData(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = "GroupVar")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "Replicate")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_id")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
rownames(metadata)<-metadata$sample_id
pheatmap(rld_cor, annotation = metadata[, c("GroupVar"), drop=F])
#rownames(metadata)<-paste0(metadata$GroupVar,"_",metadata$Replicate)
#pheatmap(rld_cor, annotation = rownames(metadata))


# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, 
               name = "GroupVar_DEGS1_vs_SCRAM",
               alpha = 0.05)

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res <- lfcShrink(dds, 
                 coef = "GroupVar_DEGS1_vs_SCRAM",
                 res=res,
                 type = "apeglm")


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

##Force Scientific notation in some columns
# res_tbl$pvalue<-formatC(res_tbl$pvalue, format = "e")
# res_tbl$padj<-formatC(res_tbl$padj, format = "e") 


###--------------        Annotate results      -----------------
## -------------------------------------------------------------

#results_annotated<-res_tbl
rownames(res_tbl)<-res_tbl$gene
# Add gene full name
res_tbl$description <- mapIds(x = org.Mm.eg.db,
                                        keys = rownames(res_tbl),
                                        column = "GENENAME",
                                        keytype = "SYMBOL",
                                        multiVals = "first")

# Add gene symbol
#res_tbl$symbol <- res_tbl$gene
rownames(res_tbl)<-res_tbl$gene
# Add ENTREZ ID
res_tbl$entrez <- mapIds(x = org.Mm.eg.db,
                                   keys = rownames(res_tbl),
                                   column = "ENTREZID",
                                   keytype = "SYMBOL",
                                   multiVals = "first")

# Add ENSEMBL
rownames(res_tbl)<-res_tbl$gene
res_tbl$ensembl <- mapIds(x = org.Mm.eg.db,
                                    keys = row.names(res_tbl),
                                    column = "ENSEMBL",
                                    keytype = "SYMBOL",
                                    multiVals = "first")
#Remove any genes that do not have any entrez identifiers

res_tbl <- subset(res_tbl, is.na(entrez) == FALSE)

res_tbl.na<-res_tbl[is.na(res_tbl$log2FoldChange),"gene"]
na_list<-res_tbl[is.na(res_tbl$log2FoldChange),"gene"]
normalized_counts.na<-normalized_counts[which(rownames(normalized_counts) %in% na_list$gene),]

## Remove genes with zero=NA form results

res_tbl <- subset(res_tbl, is.na(log2FoldChange) == FALSE)


# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

# Set thresholds
log2fc_cutoff <- 0.58

# Count significantly up/down genes above threshold
n_sig_up <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff) %>% 
  nrow()
n_sig_dn <- dplyr::filter(sig_res, log2FoldChange <= -log2fc_cutoff) %>% 
  nrow()

# Scatterplot

## Extract normalized counts from dds object
normalized_counts <- counts(dds, normalized = TRUE)

## Extract top 20 DEG from resLFC (make sure to order by padj)
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n = 42)

## Extract matching normalized count values from matrix
top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
top20_sig_counts

## Convert wide matrix to long data frame for ggplot2
top20_sig_df <- data.frame(top20_sig_counts)
top20_sig_df$gene <- rownames(top20_sig_counts)

top20_sig_df <- melt(setDT(top20_sig_df), 
                     id.vars = c("gene"),
                     variable.name = "sample_id") %>% 
  data.frame()

## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
#top20_sig_df$sample_id <- gsub("\\.", " ", top20_sig_df$sample_id)
top20_sig_df

## Join counts data frame with metadata
top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                           by = "sample_id")
top20_sig_df

## Generate plot
ggplot(top20_sig_df, aes(y = value, x = GroupVar, col = GroupVar)) +
  geom_jitter(height = 0, width = 0.15, size=1) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("condition") +
  ggtitle("Top 20 Significant DE Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene)



### Individual gene expression for Degs1 and Degs2
gene_interest<-c("Degs1","Degs2","Cers2","Smpd3","Ormdl3","Clec4g","Fabp4","Irs2","Lipg","Timp1")

#Kupffer
markers<-c("CD14,CD68,TNF,LYZ,MPO,VSIG4,PROK2,STARD5,MSR1,PPARA,TREM1,SLC11A1,MMP13,VDR,CHIT1,OSM,PPARD,CLEC4F,ADGRE1,IL1B,SLC40A1,G6PD,TIMD4,MNDA,SLC15A3,DNASE1L3,MARCO,HFE,CD38,CD163,CCR5,CLEC1B,CLEC4G,GPIHBP1,FOLR2,PLTP,FTL,IRF7,SPIC,CSF1R,C1QC,C1QA,C1QB,CLEC4E")
#Hepatocyte
markers<-c("AFP,HNF4A,KRT8,ALB,KRT18,FOSL1,EPPK1,UCP2,GCK,LRP5,SLC10A1,NOS2,ATP7B,GJB2,FGFR4,PROX1,CRP,SLC2A2,TFR2,KIF13B,LIPC,VDR,ASGR1,ARG1,G6PC,OTC,SERPINA1,ZHX2,HHEX,FOXA3,FOXA2,FOXA1,CYP2E1,CYP1A2,CEBPA,CDH1,GPC3,AHR,CPS1,GLS2,PCK1,TAT,WT1,PRRG4,SULT1A1,APOH,CTNNB1,ABCC3,FGB,AQP3,PLSCR1,FGA,APOB,ANG,ANXA13,SAT2,SFRP5,A1CF,APOA1,BNIP3,FGL1,PAH,SERPINA6,APOA2,AZGP1,FGG,APOC3,DEFB1,TM4SF4,GC,AMBP,ORM1,TTR,HAL,ASS1,SCD,HMGCS1,ACSS2,TM7SF2,SEC16B,SLBP,RND3,BCHE,GHR,ALDH6A1,MASP2,AKR1C1,HAMP,GLUL,ACLY,ASL,TMEM97,CP,SLPI,ACAT2,TM4SF5,MSMO1,LEPR,RCAN1,AR,RPP25L,HSD11B1,APOM,TKFC,G0S2,PON3,C1orf53,TTC36,FST,MCC,AQP9,GSTA2,NNT,SAA4,MRPS18C,OCIAD1,APOA5,ENTPD5,C4B,EID2,TP53INP2,ATIC,SERPINH1,SAMD5,GRB14,CD3G,RHOB,EPB41L4B,GPAT4,SPTBN1,SDC2,PHLDA1,WTAP,ACADM")
gene_interest<-str_to_title(tolower(unlist(str_split(markers,"\\,"))))

#gene_interest<-c("Degs1","Degs2")

gene_inte_tbl<-res_tbl[which(res_tbl$gene %in% gene_interest ),]

# ## Extract top 20 DEG from resLFC (make sure to order by padj)
# top_sig_genes <- gene_inte_tbl %>%
#   dplyr::arrange(padj) %>%
#   dplyr::pull(gene) %>%
#   head(n = 42)

## Extract matching normalized count values from matrix
top_sig_counts <- normalized_counts[rownames(normalized_counts) %in% gene_interest, ]
top_sig_counts

## Convert wide matrix to long data frame for ggplot2
top_sig_df <- data.frame(top_sig_counts)
top_sig_df$gene <- rownames(top_sig_counts)

top_sig_df <- melt(setDT(top_sig_df), 
                   id.vars = c("gene"),
                   variable.name = "sample_id") %>% 
  data.frame()

## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
#top_sig_df$sample_id <- gsub("\\.", " ", top_sig_df$sample_id)
top_sig_df

## Join counts data frame with metadata
top_sig_df <- plyr::join(top_sig_df, as.data.frame(colData(dds)),
                         by = "sample_id")

#top_sig_df.subset<-top_sig_df[which(top_sig_df$Time %in% filter_by),]
## Generate plot
ggplot(top_sig_df, aes(y = value, x = GroupVar, col = GroupVar)) +
  geom_jitter(height = 0, width = 0.15, size=2) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("condition") +
  ggtitle("Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +facet_wrap(~ gene, 
                                                            #ncol = 1, scales = "free_x", 
                                                            #labeller = as_labeller(
                                                            # c("Chow" = "Chow", "15weeks" = "HFD 15w", "30weeks" = "HFD 30w")), 
                                                            strip.position = "top"
                                                            # +scale_colour_manual(labels=paste0(umap_paleta$celtype," (",umap_paleta$cellalias,")"), values=umap_paleta$paleta)+
  )+ theme(
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_text(colour = "black", size = rel(0.5)),
    legend.key.size = unit(0.7, "lines"),
    legend.text = element_text(size = rel(0.7), colour = "black"),
    legend.key = element_rect(colour = NA, fill = NA),
    legend.background = element_rect(colour = NA, fill = NA),
    legend.box = "horizontal",
    panel.grid = element_blank(),
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank(),
    panel.background = element_blank(),
    #axis.line = element_blank(),
    strip.text = element_text(face = "plain")
  )+ggtitle("")+guides(colour=guide_legend(ncol=1, byrow=TRUE)) #check this it is not working

### Preparing for Gene Set Enrichment Analysis
# normalized_counts_log2<-log2(normalized_counts)
# vsd <- vst(dds, blind=FALSE)
# normalized_counts_vsd <- assay(vsd)
# rld <- rlog(dds, blind=FALSE)
# normalized_counts_rld <- assay(rld)
