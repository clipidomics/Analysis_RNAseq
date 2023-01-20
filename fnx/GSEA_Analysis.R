
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
library(openxlsx)
library(grid)
library(gridExtra)
library(svglite)
library(ggpubr)
library(genekitr)
#library(patchwork)

#Check MSigDB version
#packageVersion("msigdbr")

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
# all_gene_sets = msigdbr(species = "Mus musculus")
# rm(all_gene_sets)
# #head(all_gene_sets)
# unique(all_gene_sets$gs_cat)
# unique(all_gene_sets$gs_name)
# #rm(all_gene_sets)
# print(n=100,msigdbr_collections())

### Select Molecular Signature to Explore
MolSignatures<- c("H","C1","C2","C3","C4","C5","C6","C7","C8")
names(MolSignatures)<-c("Hallmark",
                        "Positional",
                        "Curated gen sets",
                        "Regulatory",
                        "Computational",
                        "Ontologies",
                        "Oncogenic",
                        "Immunologic",
                        "Cell type sets"
                        )
MolSignatures <-selectOption(MolSignatures,"Select Molecular Signatures to expore")
#M <- msigdbr(species = "Mus musculus", category = "C2", subcategory = c("CP:KEGG"))
Mo <-bind_rows(lapply(MolSignatures,FUN=function(x){msigdbr(species = "Mus musculus", category = x)}),.id = "column_label")[,-c(1)] 
Mo$gs_subcat[Mo$gs_subcat==""]<-"no subcat"

# Select which Subcats
MolSubcats<-unique(Mo$gs_subcat)
MolSubcats <-selectOption(MolSubcats,"Select Molecular Sub-Categories (TO REMOVE)",toRemove = T)
Mo <-Mo[which(Mo$gs_subcat %in% MolSubcats),]

unique(Mo[,c(1,2)])

## Filter specific description terms

# repeat{
#   statements...
#   if(condition){
#     break
#   }
# }

list_MolTerms<- c("List of Terms to filter separated by ,"="sphingo,ceramide,collagen,stellate,macrophage,monocyte,lipase,lipoprotein,cholesterol")
fix(list_MolTerms)

list_MolTerms<-str_split_1(list_MolTerms,pattern=",")

MolTerms<-paste0(".*(",
          paste0(paste(paste0(toupper(gsub("(.).*","[\\1|",list_MolTerms)),gsub("(.).*","\\1]",list_MolTerms),gsub("(.)(.*)","\\2",list_MolTerms)),collapse="|"),
          #paste(toupper(paste0(toupper(gsub("(.).*","[\\1|",list_MolTerms)),gsub("(.).*","\\1]",list_MolTerms),gsub("(.)(.*)","\\2",list_MolTerms))),collapse="|")
          #),
          ")"))

MolTerms
#MolTerms<-".*([S|s]phingo|[C|c]eramide|[M|m]acrophage|[M|m]onocyte|[S|s]tellate|[C|c]ollagen)|[L|l]ipase|[L|l]ipoprotein|[C|c]holesterol)"
#names(MolTerms)<-"Specific name terms to filter by: * (if all)"
#fix(MolTerms)

M<-Mo[grepl(MolTerms,Mo$gs_description),]
M_pathways_selected<-M %>% dplyr:: mutate(pathway=paste(gs_subcat,gs_name,sep="_"))%>% select(pathway,gs_description) %>% distinct()
print(M_pathways_selected,n=50)
print(paste0("You selected a list of N = ",length(M_pathways_selected$pathway)," Pathways and n_genes = ", length(unique(M$gene_symbol))))

### Load specific ontology ToxPat
Add_pathway<-"No/Yes"
names(Add_pathway)<-"Do you wish to add ToxPanel?"
fix(Add_pathway)
if (Add_pathway=="Yes"){
library(babelgene)
# CustomGeneSet<-data.table::fread(paste0(directory,"/inputs_ToxPanel_analysis/Genes_modules_associted_with_liver_damage.csv"),sep=",", header = T, skip = 0,blank.lines.skip=T)
# CustomGeneSet<-CustomGeneSet[-1,]
CustomGeneSet<-read.xlsx(xlsxFile = paste0(directory,"/inputs_ToxPanel_analysis/Genes_modules_associated_endpoints_tggates_manuscript_20150928.xlsx") , sheet = 1,sep.names = " ", skipEmptyRows = FALSE)
## Transform to list 
CustomGeneSet.ls<-as.list(CustomGeneSet)
## Purge empty genes
for (i in 1:length(CustomGeneSet.ls)){
  CustomGeneSet.ls[[i]]<-unique(CustomGeneSet.ls[[i]][CustomGeneSet.ls[[i]]!=""])
}
### Homologene Rat vs Mouse ###
# 1- Map Rat genes to Mouse. Look for orthologs between rat and mouse: Take directly symbols as thet were mouse !!!
CustomGeneSet.ls.mouse<-list()
for (i in 1:length(CustomGeneSet.ls)){
  CustomGeneSet.ls.mouse[[i]]<-orthologs(genes = CustomGeneSet.ls[[i]], species = "mouse", human = FALSE)
}
names(CustomGeneSet.ls.mouse)<-names(CustomGeneSet.ls)

## Get table
CustomGeneSet_tbl<-bind_rows(CustomGeneSet.ls.mouse, .id = "column_label")

colnames(CustomGeneSet_tbl)<-c("gs_name",
                               "human_gene_symbol",
                               "human_entrez_gene",
                               "human_ensembl_gene",
                               "taxon_id",
                               "gene_symbol",      
                               "entrez_gene", 
                               "ensembl_gene",
                               "ortholog_sources",
                               "num_ortholog_sources")

CustomGeneSet_tbl$gs_cat<-"ToxPanel"
CustomGeneSet_tbl$gs_subcat<-gsub("(.*)\\((.*)\\)","\\2",CustomGeneSet_tbl$gs_name)
CustomGeneSet_tbl$gs_description<-paste0("ToxPanel",CustomGeneSet_tbl$gs_name)
CustomGeneSet_tbl$gs_url<-"https://toxpanel.bhsai.org/toxpanel/abouttoxpanel.xhtml"

## Fusion with existing ontology TODO: if not don't fusion
Mo<-full_join(Mo,CustomGeneSet_tbl) ### to handle latter
M<-full_join(M,CustomGeneSet_tbl)
#M_pathways_selected<-unique(M[,c("gs_cat","gs_subcat","gs_name","gs_description","gs_url")])

}

#M<-CustomGeneSet_tbl
M_pathways_selected<-M %>% dplyr:: mutate(pathway=paste(gs_subcat,gs_name,sep="_"))%>% select(pathway) %>% distinct()
#print(M_pathways_selected,n=50)
print(paste0("You selected a list of N = ",length(M_pathways_selected$pathway)," Pathways and n_genes = ", length(unique(M$gene_symbol))))

#unique(M$gs_cat)
# M<-M[which(M$gs_cat %in% c("Liver")),]
##------------# Begin with GSEA analysis #-------------------########

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

padj_cutoff <- c("Fix Fold_Change padj_cutoff for Hypegeometric GSEA:"=0.05)
fix(padj_cutoff)

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

set_threshold <- c("Fix set_threshold for Hypegeometric GSEA:"= 0.2)
fix(set_threshold)


plot_Hiper_GSEA<-enrich.M.df %>% 
  filter(p.adjust <= set_threshold) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = str_sub(gsub("HALLMARK_|WP|REACTOME|KEEG","", Description),1,50),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic(10) +
  theme(axis.text.y = element_text(size=6, face="plain"))+
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = paste0("Hypergeometric enrichent of significant genes (FDR < 0.05)\nin gene sets (FDR < ",set_threshold,")")
       )

plot_Hiper_GSEA

select_type<-paste(MolSignatures, collapse="_")

Searched_pathways<-paste0("----MsigDB----\n",paste(MolSignatures,"->",MolSubcats,collapse="\n"),"\n----Filters----\n",paste(list_MolTerms,"",collapse=" \n"))

text.p <- ggparagraph(text = Searched_pathways, face = "italic", size = 11, color = "black")

ggsave(filename = paste0(directory,"/Figures/",select_type,"_",deparse(substitute(plot_Hiper_GSEA)),"_",gsub("[^A-z0-9]","",Sys.time()),".png"),
       annotate_figure(
         ggarrange(plot_Hiper_GSEA,text.p,
                   #labels = c("A", "B", "C","D","E","F"),
                   common.legend = TRUE,
                   ncol = 2, nrow = 1) #,
         #top = text_grob(paste0("Hypergeom. GSEA\n","p.adj=",padj_cutoff,"\nThreshold=",set_threshold), color = "black", face = "plain", size = 14,just="left"),
         #left=text_grob(Searched_pathways, color = "black", face = "plain", size = 6,just="left")
         ),
       width = 40, height = 20, dpi = 300, units = "cm", device='png')
browseURL(paste0(directory))
########################

#### -------------  Fold change GSEA -----------------     ####

#format gene set database
head(M)

M.ensembl.ls <- M %>% 
  dplyr::select(gs_name, ensembl_gene) %>% 
  dplyr::group_by(gs_name) %>% 
  dplyr::summarise(all.genes = list(unique(ensembl_gene))) %>% 
  deframe()

#M.ensembl.ls.copy<-M.ensembl.ls
#M.ensembl.ls<-c(CustomGeneSet.ls.ensembl,M.ensembl.ls.copy)

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
#gsea.M<- fgseaSimple(pathways = M.ensembl.ls,
                      # stats = FC.vec,
                      # scoreType = scoreType,
                      # nperm=1000)

### Check this!!!
#gsea.M<-fgsea(M.ensembl.ls, FC.vec[!duplicated(FC.vec)] , minSize=15, maxSize = 500, nperm=1000)

gsea.M<-fgsea(M.ensembl.ls, FC.vec, minSize=10, maxSize = 500, nperm=1000)

gsea.M<-gsea.M[order(gsea.M$padj),]


#gse <- genGSEA(genelist = M.ensembl.ls, geneset = gs)


### plot volcano

## volcano plot
# # get top3 of up and down pathways
# plotGSEA(gsea.M, plot_type = "volcano", show_pathway = 3)
# 
# # choose pathway by character
# pathways <- c('HALLMARK_KRAS_SIGNALING_UP','HALLMARK_P53_PATHWAY','HALLMARK_GLYCOLYSIS')
# plotGSEA(gse, plot_type = "volcano", show_pathway = pathways)


#### plot gsea results ####
#class(gsea.M)
#set_GSEA_threshold<-0.05
print(tibble(gsea.M),n=30)
set_GSEA_threshold <- c("Fix GSEA_threshold :"=0.05)
fix(set_GSEA_threshold)

pathway_formated<-gsea.M$pathway

## TODO personalize 
plot_GSEA_NES<-gsea.M %>% 
  filter(padj <= set_GSEA_threshold) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(pathway=paste0(str_to_title(sub("(.*?)\\s(.*)","\\1",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",pathway))))),
                " ",sub("(.*?)\\s(.*)","\\2",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",pathway)))))
  ) %>% 
  
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y=NES)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.y = element_text(size=16, face="plain"),
        axis.text.x = element_text(size=16, face="plain"),
        plot.title = element_text(hjust = 0.5, size = 16, face="bold"),
        axis.title.y = element_text(size = 16, margin = margin(
          t = 0,
          r = 20,
          b = 0,
          l = 0
        )),
        axis.title.x = element_text(size = 16, margin = margin(
          t = 0,
          r = 20,
          b = 0,
          l = 0
        )),
        )+
  #Force equal max min
  lims(y=c(-3.2,3.2)) +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Normalized enrichment score (NES)",
       x="Gene set",
       title = paste0("GSEA (FDR < ",set_GSEA_threshold,")\nDown <--         --> Up "))

ggsave(plot_GSEA_NES,filename = paste0(directory,"/Figures/",select_type,"_",deparse(substitute(plot_GSEA_NES)),"_",gsub("[^A-z0-9]","",Sys.time()),".svg"),
       # annotate_figure(
       #   ggarrange(plot_GSEA_NES,
       #             #labels = c("A", "B", "C","D","E","F"),
       #             common.legend = TRUE,
       #             ncol = 2, nrow = 1) #,
       #   #top = text_grob(paste0("Hypergeom. GSEA\n","p.adj=",padj_cutoff,"\nThreshold=",set_threshold), color = "black", face = "plain", size = 14,just="left"),
       #   #left=text_grob(Searched_pathways, color = "black", face = "plain", size = 6,just="left")
       # ),
       width = 30, height = 15, dpi = 300, units = "cm",device = grDevices::svg)


#svglite(filename = paste0(directory,"/Figures/",select_type,"_",deparse(substitute(plot_GSEA_NES)),"_",gsub("[^A-z0-9]","",Sys.time()),".svg"))
#plot(plot_GSEA_NES)
#dev.off()
###---------- Plot list of implicated pathways ------------

# topPathwaysUp <- gsea.M[ES > 0][head(order(pval), n=10), pathway]
# topPathwaysDown <- gsea.M[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# plotGseaTable(M.ensembl.ls[topPathways], FC.vec, gsea.M, 
#               gseaParam=0.5)

###---------- Collapse Pathways ------------ NOT WORKING !!!


collapsedPathways <- collapsePathways(gsea.M[order(pval)][padj <= set_GSEA_threshold], 
                                      M.ensembl.ls, FC.vec)
mainPathways <- gsea.M[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### Tidy names and Plot table 
list_pathways<-M.ensembl.ls[mainPathways]
names(list_pathways)<-paste0(str_to_title(sub("(.*?)\\s(.*)","\\1",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",names(list_pathways)))))),
                             " ",sub("(.*?)\\s(.*)","\\2",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",names(list_pathways))))))
gsea.M.formated<-gsea.M %>% 
    mutate(pathway=paste0(str_to_title(sub("(.*?)\\s(.*)","\\1",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",pathway))))),
                          " ",sub("(.*?)\\s(.*)","\\2",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",pathway)))))
    )

plot_GSEA_tbl<-plotGseaTable(list_pathways, FC.vec, gsea.M.formated, 
              gseaParam=0.5,
              pathwayLabelStyle=list(size=8, color="black"),
              headerLabelStyle=list(size=8, color="black"),
              valueStyle = list(size=8, color="black")
              )
plot_GSEA_tbl

### TODO: Put some plots together
#

# plot_GSEA_tbl+plot_GSEA_NES
select_type<-paste(MolSignatures, collapse="_")

#ggsave(filename = paste0(directory,"/Figures/",select_type,"_",deparse(substitute(plot_GSEA_tbl)),deparse(substitute(plot_GSEA_NES)),"_",gsub("[^A-z0-9]","",Sys.time()),".png"),
       
# ggsave(filename = paste0(directory,"/Figures/",select_type,"_",gsub("[^A-z0-9]","",Sys.time()),".png"),
#        annotate_figure(
#          ggarrange(plot_GSEA_NES,plot_GSEA_tbl,plot_Hiper_GSEA,
#                    #labels = c("A", "B", "C","D","E","F"),
#                    common.legend = TRUE,
#                    ncol = 2, nrow = 2),
#          top = text_grob(paste0("Fold Change Gene Enrichment Analysis"), color = "black", face = "plain", size = 14,just="center")
#        ),
#        width = 40, height = 25, dpi = 600, units = "cm", device='png')
ggsave(filename = paste0(directory,"/Figures/",select_type,"_",deparse(substitute(plot_GSEA_tbl)),deparse(substitute(plot_GSEA_NES)),"_",gsub("[^A-z0-9]","",Sys.time()),".svg"),
       
        annotate_figure(
         ggarrange(plot_GSEA_NES,plot_GSEA_tbl,
                   #labels = c("A", "B", "C","D","E","F"),
                   common.legend = TRUE,
                   ncol = 2, nrow = 1),
         top = text_grob(paste0("Fold Change Gene Enrichment Analysis"), color = "black", face = "plain", size = 14,just="center")
       ),
       width = 40, height = 10, dpi = 600, units = "cm", device='svg')

saveTable <- gridFtable(M_pathways_selected)
png(file = paste0(directory,"/Figures/",select_type,"_Choosen_Pathways_",gsub("[^A-z0-9]","",Sys.time()),".png"), width = saveTable$width, height = saveTable$heigth, units = "in", res = 100)
grid.newpage()
grid.table(saveTable$data, rows = NULL, theme = saveTable$theme)
dev.off()

## Write Excell sheet
write.xlsx(M_pathways_selected, file = paste0(directory,"/Figures/",select_type,"_Choosen_Pathways_",gsub("[^A-z0-9]","",Sys.time()),".xlsx"), colNames = TRUE, borders = "columns")

browseURL(paste0(directory))


#----------------CALL  HEATMAP -------------------------#########
source("fnx/GSEA_heatmap.R")
#----------------CALL  HEATMAP -------------------------#########

# ### Analysis for ToxPanel https://toxpanel.bhsai.org/toxpanel/abouttoxpanel.xhtml ###
# library(babelgene)
# # CustomGeneSet<-data.table::fread(paste0(directory,"/inputs_ToxPanel_analysis/Genes_modules_associted_with_liver_damage.csv"),sep=",", header = T, skip = 0,blank.lines.skip=T)
# # CustomGeneSet<-CustomGeneSet[-1,]
# CustomGeneSet<-read.xlsx(xlsxFile = paste0(directory,"/inputs_ToxPanel_analysis/Genes_modules_associated_endpoints_tggates_manuscript_20150928.xlsx") , sheet = 1,sep.names = " ", skipEmptyRows = FALSE)
# ## Transform to list 
# CustomGeneSet.ls<-as.list(CustomGeneSet)
# ## Purge empty genes
# for (i in 1:length(CustomGeneSet.ls)){
#   CustomGeneSet.ls[[i]]<-unique(CustomGeneSet.ls[[i]][CustomGeneSet.ls[[i]]!=""])
# }
# ### Homologene Rat vs Mouse ###
# # 1- Map Rat genes to Mouse. Look for orthologs between rat and mouse: Take directly symbols as thet were mouse !!!
# CustomGeneSet.ls.mouse<-list()
# for (i in 1:length(CustomGeneSet.ls)){
#   CustomGeneSet.ls.mouse[[i]]<-orthologs(genes = CustomGeneSet.ls[[i]], species = "mouse", human = FALSE)
# }
# names(CustomGeneSet.ls.mouse)<-names(CustomGeneSet.ls)
# 
# ## Get table
# CustomGeneSet_tbl<-bind_rows(CustomGeneSet.ls.mouse, .id = "column_label")
# 
# colnames(CustomGeneSet_tbl)<-c("gs_name",
#                                "human_gene_symbol",
#                                "human_entrez_gene",
#                                "human_ensembl_gene",
#                                "taxon_id",
#                                "gene_symbol",      
#                                "entrez_gene", 
#                                "ensembl_gene",
#                                "ortholog_sources",
#                                "num_ortholog_sources")
# 
# CustomGeneSet_tbl$gs_cat<-"ToxPanel"
# CustomGeneSet_tbl$gs_subcat<-gsub("(.*)\\((.*)\\)","\\2",CustomGeneSet_tbl$gs_name)
# CustomGeneSet_tbl$gs_description<-paste0("ToxPanel",CustomGeneSet_tbl$gs_name)
# CustomGeneSet_tbl$gs_url<-"https://toxpanel.bhsai.org/toxpanel/abouttoxpanel.xhtml"
# 
# ### 
# CustomGeneSet.ls.entrez<-list()
# for (i in 1:length(CustomGeneSet.ls.mouse)){
#   CustomGeneSet.ls.entrez[[i]]<-CustomGeneSet.ls.mouse[[i]]$entrez
# }
# names(CustomGeneSet.ls.entrez)<-names(CustomGeneSet.ls.mouse)
# ### 
# CustomGeneSet.ls.ensembl<-list()
# for (i in 1:length(CustomGeneSet.ls.mouse)){
#   CustomGeneSet.ls.ensembl[[i]]<-CustomGeneSet.ls.mouse[[i]]$ensembl
# }
# names(CustomGeneSet.ls.ensembl)<-names(CustomGeneSet.ls.mouse)
# 
# ### Transform and write result to gmt file ### 
# 
# sets_ls<-CustomGeneSet.ls.entrez
# nPaths<-length(sets_ls)
# description<-"https://toxpanel.bhsai.org/toxpanel/abouttoxpanel.xhtml"
# ###  Write a list in .gmt form  ###
# out_ls <- lapply(seq_len(nPaths), function(i){
#   c(names(sets_ls[i]),description, sets_ls[[i]])
# })
# out_char <- vapply(out_ls, paste, collapse = "\t", FUN.VALUE = character(1))
# 
# writeLines(out_char, con = paste0(directory,"/inputs_ToxPanel_analysis/Genes_modules_associated_with_liver_damage_MOUSE.gmt"))
# 
# ### Get other values for ToxPanel
# 
# ToxPanel.rts<-model.results[,c("gene","log2FoldChange","pvalue","entrez")]
# 
# ToxPanel.data<-ToxPanel.rts[c(4,2,3)]
# colnames(ToxPanel.data)<-c("Gene","logFC","pValue")
# 
# ToxPanel.data<-ToxPanel.data[!is.na(ToxPanel.data$logFC),]
#   
# write.table(ToxPanel.data, file = paste0(directory,"/inputs_ToxPanel_analysis/ToxPanel_data.txt"), 
#             quote = F,
#             sep = "\t", 
#             col.names = T,
#             row.names = F,
#             qmethod = "double")
# 
# #### Getback results from ToxPanel 
# ToxPanel.matrix<- read.xlsx(xlsxFile = paste0(directory,"/output_ToxPanel_analysis/ToxPanel_data_txt_Custom.matrix.xlsx") , sheet = 1, skipEmptyRows = FALSE)
# ToxPanel.gsea<- read.xlsx(xlsxFile = paste0(directory,"/output_ToxPanel_analysis/ToxPanel_data_txt_Custom.xlsx") , sheet = 1, skipEmptyRows = FALSE)
# 
# 
# #### PLOT heatmaps for speficific genes in ontologies
# 
# 
# 
# 
# 
# a### TODO RADAR PLOT
# 
# # install.packages("devtools")
# # devtools::install_github("ricardo-bion/ggradar")
# library(ggradar)
# df.radar<-as.numeric(t(ToxPanel.gsea[,c(7)]))
# names(df.radar)<-ToxPanel.gsea$GeneSet
# df.radar<-data.frame(df.radar)
# #df.radar$aafc.z<-as.numeric(df.radar$aafc.z)
# 
# ggradar(df.radar)
# 
# set.seed(4)
# df <- data.frame(matrix(runif(30), ncol = 10))
# df[, 1] <- paste0("G", 1:3)
# colnames(df) <- c("Group", paste("Var", 1:9))
# ggradar(df)
# 
# 
# ### Running GSEA for ToxPanel pathways
# 
# 
