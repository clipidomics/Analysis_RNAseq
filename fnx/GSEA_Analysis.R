# #### Install packages ####
# install.packages("msigdbr")
# install.packages("fgsea")
# 
# install.packages("BiocManager")
# BiocManager::install("clusterProfiler")

#### Load packages ####
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)

#Check MSigDB version
packageVersion("msigdbr")

#### Load data ####
#load("data/RSTR_data_clean_subset.RData")

#model.results <- read_csv("data/RSTR.Mtb.model.subset.csv")

model.results <- res_tbl

### Check some duplicated genes, genes not annotated, without ensembl

ifelse(sum(duplicated(as.character(model.results$gene)))>0,"CHECK FOR DUPLICATED GENES !!","GREAT, NOT DUPLICATED GENES !!")
ifelse(sum(duplicated(as.character(model.results$ensembl)))>0,"CHECK FOR DUPLICATED Ensmbl GENES !!","GREAT, NOT DUPLICATED GENES !!")

#look for genes without ensembl, Check lattest version of org.Mm.eg.db
model.results.NA<-model.results[is.na(model.results$ensembl),]
# genes starting with 'Gm-' or ending with 'Rik' include protein-coding RNAs, long non-coding RNAs, and antisense transcripts.
model.results.NA_not_RikGm<-model.results.NA[!grepl("Rik|Gm",model.results.NA$gene),]
#length(model.results.NA_not_RikGm$gene)
## Number of genes Rik or Gm
model.results.Rik_Gm<-model.results[grepl("Rik|Gm",model.results$gene),]
print(paste0("Number of genes without annotation in ensembl = ",length(model.results.NA$gene)))
print(paste0("Not annotated and not .*Rik or Gm.* = ",length(model.results.NA_not_RikGm$gene)))
print(paste0("total.*Rik or Gm.* in model.results = ",length(model.results.Rik_Gm$gene)))

#### --------- Enrichment -------------- ####

#Get gene set database

msigdbr_species()
#all_gene_sets = msigdbr(species = "Mus musculus")
#head(all_gene_sets)
#unique(all_gene_sets$gs_cat)
#rm(all_gene_sets)
print(n=100,msigdbr_collections())

#M <- msigdbr(species = "Mus musculus", category = "C2", subcategory = c("CP:KEGG"))
Mo <- msigdbr(species = "Mus musculus", category = "C2")
unique(Mo$gs_subcat)
# unique(M$gs_description)

## Select specific subcategories
M<-Mo[grepl(".*(REACTOME|WIKIPATHWAYS)",Mo$gs_subcat),]

## Select specific description terms
M<-Mo[grepl(".*([S|s]phingo|[C|c]eramide)",Mo$gs_description),]
M<-Mo[grepl(".*([S|s]phingo|[C|c]eramide|[M|m]acrophage|[M|m]onocyte|[S|s]tellate|[C|c]ollagen)|[L|l]ipase|[L|l]ipoprotein|[C|c]holesterol|[G|g]augher",Mo$gs_description),]

class(M)

#Define significant genes
##Look at FDR spread
ggplot(model.results, aes(x=padj)) +
  geom_histogram(bins=100) +
  theme_classic()

table(model.results$padj == 0)

# #Set FDR cutoff and get genes 
# signif <- model.results %>% 
#   filter(FDR <= 1E-16)


# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
signif <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)


#get entrez ID
signif.entrez <- unique(signif$entrez)
M.entrez <- select(M, gs_name, entrez_gene)

# -----------     Run hypergemometric enrichment ----------------

enrich.M <- enricher(gene = signif.entrez, TERM2GENE = M.entrez)

# #or for ensembl
# #get ensembl ID
# signif.ensembl <- unique(signif$ensembl_gene_id)
# H.ensembl <- select(H, gs_name, ensembl_gene)
# 
# #run enrichment
# enrich.H <- enricher(gene = signif.ensembl, TERM2GENE = H.ensembl)

#### Extract results ####
class(enrich.M)

head(enrich.M@result)

class(enrich.M@result$GeneRatio)

#format results
enrich.M.df <- enrich.M@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)

#### Visualize results ####
set_threshold<-0.05
enrich.M.df %>% 
  filter(p.adjust <= set_threshold) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = gsub("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = paste0("Hypergeometric enrichent of significant genes (FDR < 0.05)\nin gene sets (FDR < ",set_threshold,")")
       )


########################

#### -------------  Fold change GSEA -----------------     ####

#format gene set database
head(M)

M.ensembl.ls <- M %>% 
  select(gs_name, ensembl_gene) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(ensembl_gene))) %>% 
  deframe()

####  Calculate fold change ### 
### it is not necessary DESeq2 already calculate this parameter
#Extract expression data
# FC <- as.data.frame(dat$E) %>% 
#   #Move gene IDs from rownames to a column
#   rownames_to_column("ensembl_gene_id") %>% 
#   #Make long format with all expression values in 1 column
#   pivot_longer(-ensembl_gene_id, 
#                names_to = "libID", values_to = "expression") %>% 
#   #Extract RSID and TB condition from libID
#   #If this info was not in the libID, we could get it by joining
#   # with dat$targets
#   separate(libID, into = c("RSID","condition"), sep="_") %>% 
#   #Make wide with media and tb expression in separate columns
#   pivot_wider(names_from = condition, values_from = expression) %>% 
#   #Calculate tb minus media fold change (delta for change)
#   #Because expression values are log2, subtraction is the same as division
#   mutate(delta = TB-MEDIA) %>% 
#   #Calculate mean fold change per gene
#   group_by(ensembl_gene_id) %>% 
#   summarise(mean.delta = mean(delta, na.rm=TRUE)) %>% 
#   #Arrange by descending fold change
#   arrange(desc(mean.delta))

##format for gsea 

# FC.vec <- FC$mean.delta
# names(FC.vec) <- FC$ensembl_gene_id

FC.vec <- model.results$log2FoldChange
names(FC.vec) <- model.results$ensembl

#Check no NAs in log2FoldChange

if(sum(is.na(names(FC.vec)))>0){
  print(paste0("Warning: FC.vec ranked vector contain ",sum(is.na(names(FC.vec))), " non-annotated genes, they will be excluded from GSEA analysis"))
  FC.vec<-FC.vec[!is.na(names(FC.vec))]
}
if(sum(is.na(FC.vec))>0){
  print(paste0("Warning: FC.vec ranked vector contain ",sum(is.na(FC.vec)), " Log2FC NAs values, they will be excluded from GSEA analysis"))
  FC.vec<-FC.vec[!is.na(FC.vec)]
}



#Check no NAs in log2FoldChange
#sum(is.na(model.results$log2FoldChange)>0)

#set score type
min(FC.vec)
max(FC.vec)
scoreType <- "std"

#Run GSEA
gsea.M<- fgseaSimple(pathways = M.ensembl.ls,
                      stats = FC.vec,
                      scoreType = scoreType,
                      nperm=1000)

### Check this!!!
#gsea.M<-fgsea(M.ensembl.ls, FC.vec[!duplicated(FC.vec)] , minSize=15, maxSize = 500, nperm=1000)

gsea.M<-fgsea(M.ensembl.ls, FC.vec, minSize=15, maxSize = 500, nperm=1000)



#### plot gsea results ####
class(gsea.H)
set_GSEA_threshold<-0.05

## TODO personalize 
plot_GSEA_NES<-gsea.M %>% 
  filter(padj <= set_GSEA_threshold) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(pathway = str_to_title(tolower(gsub("HALLMARK_|WP|REACTOME|KEEG","", pathway))),
         pathway = gsub("_"," ", pathway)) %>% 
  
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y=NES)) +
  geom_col() +
  theme_classic() +
  #Force equal max min
  lims(y=c(-3.2,3.2)) +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Normalized enrichment score (NES)",
       x="Gene set",
       title = paste0("GSEA (FDR < ",set_GSEA_threshold,")\nDown <--         --> Up "))


###---------- Plot list of implicated pathways ------------

topPathwaysUp <- gsea.M[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- gsea.M[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(M.ensembl.ls[topPathways], FC.vec, gsea.M, 
              gseaParam=0.5)

###---------- Collapse Pathways ------------ NOT WORKING !!!


collapsedPathways <- collapsePathways(gsea.M[order(pval)][padj <= set_GSEA_threshold], 
                                      M.ensembl.ls, FC.vec)
mainPathways <- gsea.M[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### Tidy names and Plot table 
list_pathways<-M.ensembl.ls[mainPathways]
names(list_pathways)<-str_to_title(tolower(gsub("HALLMARK_|WP|REACTOME|KEEG","", names(list_pathways))))
names(list_pathways)<-gsub("_"," ", names(list_pathways))

gsea.M.formated<-gsea.M %>% 
    mutate(pathway = str_to_title(tolower(gsub("HALLMARK_|WP|REACTOME|KEEG","", pathway))),
           pathway = gsub("_"," ", pathway))

plot_GSEA_tbl<-plotGseaTable(list_pathways, FC.vec, gsea.M.formated, 
              gseaParam=0.5,
              pathwayLabelStyle=list(size=8, color="black"),
              headerLabelStyle=list(size=8, color="black"),
              valueStyle = list(size=8, color="black")
              )
plot_GSEA_tbl

### TODO: Put some plots together
library(patchwork)
# library(ggpubr)
# plot_GSEA_tbl+plot_GSEA_NES
select_type<-"GSEA_Plot_combo"

ggsave(filename = paste0(directory,"/Figures/",select_type,"_",gsub("[^A-z0-9]","",Sys.time()),".png"),
       annotate_figure(
         ggarrange(plot_GSEA_NES,plot_GSEA_tbl,
                   #labels = c("A", "B", "C","D","E","F"),
                   common.legend = TRUE,
                   ncol = 1, nrow = 2),
         top = text_grob(paste0("Fold Change Gene Enrichment Analysis"), color = "black", face = "plain", size = 14,just="center")
       ),
       width = 20, height = 22, dpi = 600, units = "cm", device='png')




