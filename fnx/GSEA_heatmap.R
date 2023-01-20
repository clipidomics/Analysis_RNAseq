
normalizeTo1<-function(inputData,normalizebygroup){
  
  inputData=data.frame(inputData,stringsAsFactors = F)
  ## creamos grupos
  temp=cbind(ID=rownames(inputData),
             GROUP=gsub("(^.*)\\_(.*)","\\1",rownames(inputData)),
             inputData)
  ## seleccionamos el que vamos a utilizar para normalizar
  if (missing(normalizebygroup)){
    normalizeGroup=selectOption(as.character(unique(temp$GROUP)),name="Reference")
    normalizeGroup=paste0(normalizeGroup,collapse="|")
  } else (normalizeGroup=normalizebygroup)
  ## definimos la funcion que nos va a hacer la normalizaci?n con filas correspondientes al grupo seleccionado
  Norm_by_group=function(colValues,
                        refRow=grep(normalizeGroup,temp$ID)){
    colValues=colValues*(1/my.mean(colValues[refRow]))
    return(colValues)
  }
  #aplicamos la transformacion
  temp[,-c(1,2)]=temp%>%summarise_if(is.numeric,Norm_by_group)
 return(temp%>%select(-c("ID","GROUP")))
}


### COMPLEX HEATMAP ####

# library(devtools) 
# install_github("jokergoo/ComplexHeatmap") #instalar complexheatmap
library(ComplexHeatmap)
library(cluster)
library(circlize) #generar colores

####
# collapsedPathways <- collapsePathways(gsea.M[order(pval)][padj <= set_GSEA_threshold], 
#                                       M.ensembl.ls, FC.vec)
# mainPathways <- gsea.M[pathway %in% collapsedPathways$mainPathways][
#   order(-NES), pathway]


gsea.M.sig<-gsea.M %>% #filter(pathway %in% mainPathways)   %>% 
  filter(padj <= set_GSEA_threshold)  ##Filltro genes en pathways significativas
# Get genes associated with pathway
M.in_pathway<-M[M$gs_name %in% gsea.M.sig$pathway,] 
unique(M.in_pathway$gs_name)

# Get log2FC and other values from model.result table
M.in_pathway.join<-inner_join(M.in_pathway,model.results,by=c("gene_symbol"="gene"))
#rownames(M.in_pathway.join)<-paste0(rownames(M.in_pathway.join),"_",M.in_pathway.join$gene_symbol)

## Select 

#M.in_pathway.subset<-M.in_pathway.join[grep("SPHINGOLIPID|CERAMIDE",M.in_pathway.join$gs_name),]

#M.in_pathway.subset<-M.in_pathway.join[grep("CHOLESTEROL|SPHINGOLIPID|CERAMIDE|MONOCYTE|MACROPHAGE|INFLAMMATION|COLLAGEN|TRIGLYCERIDE",M.in_pathway.join$gs_name),]
M.in_pathway.subset<-M.in_pathway.join

#M.in_pathway.subset<-M.in_pathway.join[grep("Tox",M.in_pathway.join$gs_cat),]

# Set thresholds
FC_padj<-0.05

# Subset the significant results
M.in_pathway.subset.sig <- dplyr::filter(M.in_pathway.subset, padj < FC_padj) %>%
  dplyr::arrange(padj)
#M.in_pathway.subset.sig<-M.in_pathway.subset
# # Set thresholds for log2fc_cutoff
# Heat_log2fc_cutoff <- 0.58
# M.in_pathway.subset.sig <- dplyr::filter(M.in_pathway.subset, log2FoldChange >= Heat_log2fc_cutoff | log2FoldChange <= -Heat_log2fc_cutoff) %>%
#   dplyr::arrange(padj)

# Check significant genes output
M.in_pathway.subset.sig 
unique(M.in_pathway.subset.sig$gs_name)
## Get list of genes in pathway
Genes_sig_in_pathway<-c(M.in_pathway.subset.sig$gene_symbol)
Genes_sig_in_pathway
#subset<-data.frame(M.in_pathway.join[M.in_pathway.join$gs_cat %in% c("ToxPanel"),])

###HEATMAP
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)
mat.all<-assay(ddsMat_rlog)

mat<-mat.all[rownames(mat.all) %in% Genes_sig_in_pathway,]


#------------------------ BEGIN ANALYSIS HEATMAP ----------------------------


inputData<-data.frame(t(mat))

##PLAYING WITH TWO LEVELS

heatmap_input<-inputData

#### sumas
dfc.dots=inputData
#dfc.dots$tto = gsub("(\\W*)_.*", "\\1", rownames(dfc.dots))
dfc.dots$tto = gsub("(\\W*_\\d+)_.*", "\\1", rownames(dfc.dots))
dfc.dots = reshape2::melt(dfc.dots, id = c("tto"))
#dfc.dots$variable = gsub("(.*)_(\\d+$)", "\\1", dfc.dots$variable)
#dfc.dots = ddply(dfc.dots, .(variable, tto), colwise(myfun2, c("value")))
dfc.dots = reshape2::dcast(dfc.dots, tto ~ variable)
rownames(dfc.dots) = dfc.dots[,1]
heatmap_input = dfc.dots[,-1]


ask_specie=c("Prefijo para la imagen: "="Raton")
#fix(ask_specie)
###MATCH CASES WITH SECOND CRITERIA###

#juegas con la escala de colores
heatmap_base_color=as.data.frame(
  cbind("Parameters"=c("Heatmap color value min","Heatmap color value mean","Heatmap color value max",
                       "Heatmap colors gradient min","Heatmap colors gradient mean","Heatmap colors gradient max"),
        "Input"=c(c(-2, 0, 2),c("red", "transparent", "blue"))))
fix(heatmap_base_color)
heatmap_base_color$Input=as.character(heatmap_base_color$Input)
col_fun = colorRamp2(as.numeric(heatmap_base_color$Input[1:3]), c(heatmap_base_color$Input[4:6]))


### restart from here  

heatmap_input = dfc.dots[,-1]


  ## results id removing batch digit at the end  
  #matchTo=gsub("(.*)\\D(\\d+$)","\\2",colnames(t(heatmap_input)))
  matchTo=colnames(heatmap_input)
  
  table(matchTo)
  
  ## matching to database
  BBDD_temp=data.frame(M.in_pathway.subset.sig)
  #Filter BBDD_temp for genes in heatmap_input only
  BBDD_temp<-BBDD_temp[which(BBDD_temp$gene_symbol %in% colnames(heatmap_input)),]
  unique(BBDD_temp$gene_symbol) %in% colnames(heatmap_input)
  # BBDD_genSig<-unique(BBDD_temp$gene_symbol)
  # heatmap_genes<-colnames(heatmap_input)
  # 
  # BBDD_genSig %in% heatmap_genes
  # heatmap_genes %in%  BBDD_genSig 
  
  ## Filter heatmpa fo  genes in BBDD_temp only 
  #heatmap_input<-heatmap_input[,!is.na(match(unique(BBDD_temp$gene_symbol),colnames(heatmap_input)))]
  heatmap_input<-heatmap_input[,unique(BBDD_temp$gene_symbol)]
  
  if(sum(unique(BBDD_temp$gene_symbol) %in% colnames(heatmap_input)==colnames(heatmap_input) %in% unique(BBDD_temp$gene_symbol))>0) print("------TODO CORRECTO------") else ("REVISA heatmap_input &  M.in_pathway.subset.sig")
  
  ### FILTRAR LOS GRUPOS
  sel_columnToFilterGroups=selectOption(colnames(BBDD_temp),"ColumnForGroupFilter",toRemove=F)
  #BBDD_sel_tto=selectOption(as.character(unique(gsub("(.*)\\D(\\d+$)","\\1",BBDD_temp[,sel_columnToFilterGroups]))),name="GroupsToRemove",toRemove=T)
  BBDD_sel_tto=selectOption((unique(BBDD_temp[,sel_columnToFilterGroups])),name="GroupsToRemove",toRemove=T)
  
  if(length(BBDD_sel_tto)>0){
    BBDD_temp<-BBDD_temp[c(BBDD_temp[,sel_columnToFilterGroups]) %in% c(BBDD_sel_tto),]
  }
  
  #Format names
  BBDD_temp[,sel_columnToFilterGroups]<-paste0(str_to_title(sub("(.*?)\\s(.*)","\\1",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",BBDD_temp[,sel_columnToFilterGroups])))))," ",
  sub("(.*?)\\s(.*)","\\2",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",BBDD_temp[,sel_columnToFilterGroups])))))
  
  
  anolist<-list()
  catlist<-unique(BBDD_temp[,sel_columnToFilterGroups])
  catlist
  # ##### Ask if Collapse categories ######
  unique(BBDD_temp$gs_name)
  # ##### Collapse categories TODO ######
  ask_collapse<-data.frame(cbind(Actual=catlist),fusion="")
  #rownames(ask_collapse)<-catlist
  fix(ask_collapse)
  ask_collapse[,1]
  var_fusion<-ask_collapse[,2]
  var_fusion<-unique(var_fusion[grepl("\\S",var_fusion)])
  var_fusion
  
  ### DEBUG FROM HERE   ####
  ## Variables
  # BBDD_temp_copy<-BBDD_temp
  # mat_copy<-mat
  # heatmap_input_copy<-heatmap_input
  # heatmap_base_color_copy<-heatmap_base_color
  # sel_columnToFilterGroups_copy<-sel_columnToFilterGroups
  #   
  # BBDD_temp<-BBDD_temp_copy
  # mat<-mat_copy
  # heatmap_input<-heatmap_input_copy
  # heatmap_base_color<-heatmap_base_color_copy
  # sel_columnToFilterGroups<-sel_columnToFilterGroups_copy


  
  for (i in 1:length(var_fusion)){
    to_collapse<-ask_collapse[which(ask_collapse$fusion==var_fusion[i]),]
    BBDD_temp[BBDD_temp[,sel_columnToFilterGroups] %in% to_collapse[,1],sel_columnToFilterGroups]<-paste0(to_collapse[,1],sep="\n",collapse="")
  }
  length(unique(BBDD_temp[,sel_columnToFilterGroups]))
  #delete character at end of string
  BBDD_temp[,sel_columnToFilterGroups]<-sub("\n$","",BBDD_temp[,sel_columnToFilterGroups])
  
  #re-define catlist
  catlist<-unique(BBDD_temp[,sel_columnToFilterGroups])
  catlist
  #re-define heatmap_input
  unique(BBDD_temp$gene_symbol)
  
  #heatmap_input[,unique(BBDD_temp$gene_symbol) %in% colnames(heatmap_input)]
  heatmap_input<-heatmap_input[,unique(BBDD_temp$gene_symbol)]
  #heatmap_input<-heatmap_input[,!is.na(match(unique(BBDD_temp$gene_symbol),colnames(heatmap_input)))]
  
  #catlist<-catlist[-1]
  #dev.new()
  display.brewer.all()
  ask_annotation<-c("Elegir Set de colores:"="Paired/Dark2/Set2/Pastel2")
  fix(ask_annotation)
  ask_annotation<-unlist(str_split(ask_annotation,pattern = "/"))
  
  aux_df<-list()
  
  for(j in 1:length(catlist)){
    #j=1
    ### filtering by each posible category from selected column
    #x=unique(BBDD_temp[,sel_columnToFilterGroups])[2]
    col_df=c(sel_columnToFilterGroups,"gene_symbol")
   # df<-BBDD_temp[BBDD_temp$gene_symbol %in% colnames(heatmap_input),c(3,4)]
    df<-BBDD_temp[BBDD_temp$gene_symbol %in% colnames(heatmap_input),col_df]
    
    #distinct(df$gene_symbol)
    df<- df %>% distinct(gene_symbol, .keep_all = TRUE)
    
    ## restoring categories in case eliminated during duplicate filter
    df[paste0(df$gene_symbol,"_",catlist[j]) %in% paste0(BBDD_temp$gene_symbol,"_",BBDD_temp[,sel_columnToFilterGroups]),1]<-catlist[j]
    
    table(df[,2])
    df[df[,1]!=catlist[j],1]<-""
    # 
    # auxi<-data.frame(cbind(gene_symbol=colnames(heatmap_input),""))
    # colnames(auxi)[2]<-sel_columnToFilterGroups
    # df<-left_join(df,auxi,keep=F)
    # df<- df %>% distinct(gene_symbol, .keep_all = TRUE)
    # # 
    
    #colnames(heatmap_input) %in% df[,2]
  
    df<-df[match(colnames(heatmap_input),df$gene_symbol),] #PROBLEMATIC LINE INDUCE ERRORS!!!

    ##Auxiliar table
    #aux_df<-list()
    aux_df[[j]]<-c(df)
    names(aux_df)[j]<-catlist[j]
    #class_list=c(unique(BBDD_temp[,sel_columnToFilterGroups]))
    # Easier alternative
    aux_scale<-unlist(lapply(seq_len(length(ask_annotation)),function(i){
      c(brewer.pal(brewer.pal.info[ask_annotation[i],]$maxcolors,ask_annotation[i]))}))
    escala_x1<-c("white",aux_scale[j])
    names(escala_x1)=c(unique(df[,1])[order(unique(df[,1]))])
    #names(escala_x1)[2]<-str_to_title(tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",names(escala_x1)[2]))))
    #escala_x1[which(names(escala_x1)!=ask_annotation[2,2])]="white"
    ## saving annotation    
    anolist[[j]]<-rowAnnotation(Category = df[,1],
                          col = list(Category=escala_x1),
                          show_annotation_name = T, #False if not rows 
                          na_col = "white",
                          #name = ask_annotation[2,2],
                          name = names(escala_x1)[2],
                          #name = j,
                          #width=unit(5,"mm"),
                          #annotation_label=catlist[j],
                          annotation_label=as.character(j),
                          annotation_name_rot = 0,
                          annotation_name_side = "top"
                          )
    #anolist[[j]]@name<-str_to_title(tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",anolist[[j]]@name))))
    
    
    
    
  }

 names(aux_df)
 Cat_table<-lapply(seq_len(length(aux_df)), function(i){
   aux_df[[i]][[1]]
 })
 Cat_table<-data.frame(cbind(gene_symbol=aux_df[[1]][[2]],do.call(cbind,Cat_table)))
 
  
construct_annotation<-unlist(lapply(seq_len(length(anolist)), function(i){
    paste0("anolist[[",i,"]]")
  }))


#### HEATMAP SIMPLE #####
selection=selectOption(unique(gsub("(.*)\\D(\\d+$)","\\1",rownames(heatmap_input))),"toRemove",toRemove = T)
  
#heatmap_input = dfc.dots[,-1]
heatmap_input = t(scale((heatmap_input)))

rownames(heatmap_input) %in% aux_df[[1]][[2]] 

#Adjust number of rows to selected gene_name categories TODO: See what happened with other criteria  
#heatmap_input<-heatmap_input[rownames(heatmap_input) %in% BBDD_temp[BBDD_temp$gs_name %in% catlist,"gene_symbol"],]

#heatmap_input<-heatmap_input[rownames(heatmap_input) %in% aux_df[[1]][[2]],]
#match(rownames(heatmap_input), aux_df[[1]][[2]])

#### ordenar las clases

# asklevels=c("Cambiar order clases?"="y")
# fix(asklevels)
# 
# if(tolower(asklevels)=="y" | !exists("heatmap_levels")){
#   levels=unique(gsub("(.*)\\D(\\d+$)","\\1",rownames(heatmap_input)))
#   if(!exists("heatmap_levels")){
#     heatmap_levels=as.data.frame(cbind(levels=levels,order=c(1:length(levels))),stringAsFactor=F)
#   } else if(sum(heatmap_levels$levels %in% levels)!=length(levels)){
#     heatmap_levels=as.data.frame(cbind(levels=levels,order=c(1:length(levels))),stringAsFactor=F)
#   }
#   fix(heatmap_levels)
# }

### ordenar los grupos

asklevels_group=c("Cambiar order grupos?"="y")
fix(asklevels_group)

if(tolower(asklevels_group)=="y" | !exists("heatmap_levels_group")){
  groups=unique(gsub("(.*)_(.*)","\\1",colnames(heatmap_input)))
  if(!exists("heatmap_levels_group")){
    heatmap_levels_group=as.data.frame(cbind(levels=groups[order(groups)],order=1:length(groups)),stringAsFactor=F)
  } else if(sum(heatmap_levels_group$levels %in% groups)!=length(groups)){
    heatmap_levels_group=as.data.frame(cbind(levels=groups[order(groups)],order=1:length(groups)),stringAsFactor=F)
  }
  fix(heatmap_levels_group)
}

heatmap_input[is.na(heatmap_input)]=0
heatmap_input[is.nan(heatmap_input)]=0


plot_htmp<-Heatmap(heatmap_input,
          name = "log2FC",
          width=ncol(heatmap_input)*unit(7,"mm"),
          height=nrow(heatmap_input)*unit(3,"mm"),
          border_gp = gpar(col = "darkgrey", lty = 1,alpha=1),
          rect_gp = gpar(col = "white", lwd = .2,lty = 1),
          #column_title = "Heatmap rat?n", ## quita titulos del split
          col = col_fun,
          #ordenamiento
          
          #cluster_columns = cluster_within_group(heatmap_input, group),
          #cluster_column_slices = F,
          #cluster_row_slices = F,## tiene que ser false para poder ordenar los grupos 
          #clustering_method_rows = "median",#centroid,median,mcquitty,average,complete,single,ward.D2,ward.D
          #cluster_rows = function(m) as.dendrogram(diana(m)),
          #cluster_columns = function(m) as.dendrogram(diana(m)),
          #cluster_rows = T,
          #clustering_distance_rows = "pearson",
          #name = "mat", 
          cluster_rows = diana,
          cluster_columns=F,
          clustering_method_columns="median",
          #column_km = 3,
          #nombres filas/columnas
          show_row_dend  = F,
          column_dend_reorder = T, #controla ordenamiento del dendrograma
          show_row_names = T,
          show_column_names = F, 
          show_column_dend = F,
          # parametros split
          # row_split = factor(rownames(heatmap_input),
          #                    levels=as.character(heatmap_levels$levels)[match(rownames(heatmap_levels),heatmap_levels$order)],ordered = T),
          
          #column_split = factor(gsub("(.*)\\D(\\d+$)","\\1",colnames(heatmap_input))),
          # column_split = factor(gsub("(.*)\\D(\\d+$)","\\1",colnames(heatmap_input)),
          #                       levels=as.character(heatmap_levels_group$levels)[match(rownames(heatmap_levels_group),heatmap_levels_group$order)],ordered = T),
          column_split=factor(gsub("(.*)_(\\d+$)","\\1",colnames(heatmap_input)),
                              levels=as.character(heatmap_levels_group$levels)[match(rownames(heatmap_levels_group),heatmap_levels_group$order)],ordered = T),
          column_title_side = "top",
          #column_title_gp = gpar(fontsize = 10), #Cajitas para los titulos
          row_gap = unit(c(1), "mm"),
          column_gap = unit(c(2), "mm"),
          ### etiquetas de casos en filas
          #right_annotation  = c(anolist[[1]]),
          right_annotation  = eval(parse(text=paste0("c(",paste(construct_annotation,collapse=","),")"))), # to build annotation
          heatmap_legend_param = list(
            legend_direction = "horizontal", 
            legend_width = unit(6, "cm")),
          column_names_rot = c(45),
          column_title_rot = 0,
          column_title_gp = gpar(fontsize = 15),
          row_title_gp = gpar(fontsize = 15),
          column_names_gp = gpar(fontsize = 10),
          row_names_gp = gpar(fontsize = 10),
          row_title_rot = 0
          
  )
heigth_heatmap<-length(rownames(heatmap_input))*0.55
png(paste0(directory,"/Figures/",select_type,"_",deparse(substitute(htmp)),"_",gsub("[^A-z0-9]","",Sys.time()),".png"),
  #paste0(directory,"/Figures/HeatMap_aux_test_",as.numeric(Sys.time()),".png"),
    width = 40, height = heigth_heatmap, units = "cm", res = 200)
print(draw(plot_htmp, heatmap_legend_side="bottom", annotation_legend_side="right",
           legend_grouping = "original")
)
dev.off()

#Save as .svg file

svglite(filename=paste0(directory,"/Figures/",select_type,"_",deparse(substitute(htmp)),"_",gsub("[^A-z0-9]","",Sys.time()),".svg"),
        width = 20/2.54, height = heigth_heatmap/2.54, scaling=1)
        plot(draw(plot_htmp, heatmap_legend_side="bottom", annotation_legend_side="top",
             legend_grouping = "original")  
             )
dev.off()


# ff<-draw(plot_htmp, heatmap_legend_side="bottom", annotation_legend_side="right",
#      legend_grouping = "original")
# ggsave( ff,filename = paste0(directory,"/Figures/",select_type,"_",deparse(substitute(htmp)),"_",gsub("[^A-z0-9]","",Sys.time()),".svg"),
#        width = 40, height = 10, dpi = 600, units = "cm", device='svg')      

browseURL(paste0(directory))

#--------------CALL INDIVIDUAL GENES PLOT --------------##
source("fnx/GSEA_plot_individual_genes.R")
#--------------CALL INDIVIDUAL GENES PLOT --------------##
#rownames(Cat_table)<-aux_df[[1]][[2]]


