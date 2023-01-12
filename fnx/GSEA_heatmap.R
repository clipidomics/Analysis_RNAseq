#preparamos los datos
normalizeTo1 = function(inputData){
  
  inputData=data.frame(inputData)
  ## creamos grupos
  temp=cbind(ID=rownames(inputData),
             GROUP=gsub("(^.*)\\_(.*)","\\1",rownames(inputData)),
             inputData)
  ## seleccionamos el que vamos a utilizar para normalizar
  normalizeGroup=selectOption(as.character(unique(temp$GROUP)),name="Reference")
  normalizeGroup=paste0(normalizeGroup,collapse="|")
  ## definimos la funcion que nos va a hacer la normalizaci?n con filas correspondientes al grupo seleccionado
  normalizeTo1=function(colValues,
                        refRow=grep(normalizeGroup,temp$ID)){
    colValues=colValues*(1/my.mean(colValues[refRow]))
    return(colValues)
  }
  #aplicamos la transformacion
  temp[,-c(1,2)]=temp%>%summarise_if(is.numeric,normalizeTo1)
  return(temp[,-c(1,2)])
}


### COMPLEX HEATMAP ####

# library(devtools) 
# install_github("jokergoo/ComplexHeatmap") #instalar complexheatmap
library(ComplexHeatmap)
library(cluster)
library(circlize) #generar colores


#manual
#browseURL("https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#horizon-chart-annotation")





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
  
  ### FILTRAR LOS GRUPOS
  sel_columnToFilterGroups=selectOption(colnames(BBDD_temp),"ColumnForGroupFilter",toRemove=F)
  #BBDD_sel_tto=selectOption(as.character(unique(gsub("(.*)\\D(\\d+$)","\\1",BBDD_temp[,sel_columnToFilterGroups]))),name="GroupsToRemove",toRemove=T)
  BBDD_sel_tto=selectOption((unique(BBDD_temp[,sel_columnToFilterGroups])),name="GroupsToRemove",toRemove=T)
  
  if(length(BBDD_sel_tto)>0){
    BBDD_temp<-BBDD_temp[c(BBDD_temp[,sel_columnToFilterGroups]) %in% c(BBDD_sel_tto),]
  }
  
  anolist<-list()
  catlist<-unique(BBDD_temp[,sel_columnToFilterGroups])
  
  # ##### Collapse categories TODO ######
  # ask_collapse<-data.frame(cbind(Actual=catlist),fusion="")
  # #rownames(ask_collapse)<-catlist
  # fix(ask_collapse)
  # to_collapse<-ask_collapse[which(ask_collapse$fusion!=""),]
  # #unique(BBDD_temp[which(BBDD_temp[,sel_columnToFilterGroups] %in% to_collapse),c("gene_symbol")])
  # #BBDD_temp[which(ask_collapse$fusion!=""),sel_columnToFilterGroups]<-
  # BBDD_temp[to_collapse[,1],sel_columnToFilterGroups]<-paste0("FUSION:",BBDD_temp[which(ask_collapse$fusion!=""),sel_columnToFilterGroups],sep="\n",collapse="")
  # BBDD_temp[to_collapse[,1],sel_columnToFilterGroups]<-paste0("FUSION:",BBDD_temp[which(ask_collapse$fusion!=""),sel_columnToFilterGroups],sep="\n",collapse="")
  # 
  #catlist<-catlist[-1]
  ask_annotation<-c("Elegir Set de colores:"="Paired/Set1/Set2")
  fix(ask_annotation)
  ask_annotation<-unlist(str_split(ask_annotation,pattern = "/"))
  
  for(j in 1:length(catlist)){
    #j=1
    ### filtering by each posible category from selected column
    #x=unique(BBDD_temp[,sel_columnToFilterGroups])[2]
    df<-BBDD_temp[BBDD_temp$gene_symbol %in% colnames(heatmap_input),c(2,4)]
    #distinct(df$gene_symbol)
    df<- df %>% distinct(gene_symbol, .keep_all = TRUE)
    
    ## restoring categories in case eliminated during duplicate filter
    df[paste0(df$gene_symbol,"_",catlist[j]) %in% paste0(BBDD_temp$gene_symbol,"_",BBDD_temp[,sel_columnToFilterGroups]),1]<-catlist[j]
    
    table(df[,2])
    df[df[,1]!=catlist[j],1]<-""
    
    #
    auxi<-data.frame(cbind(gene_symbol=colnames(heatmap_input),gs_subcat=""))
   # df<-full_join(df,auxi,keep=F)
    df<-left_join(df,auxi,keep=F)
    df<- df %>% distinct(gene_symbol, .keep_all = TRUE)
    
    #df<-df[match(colnames(heatmap_input),df$gene_symbol),]
    
    #df<-df[match(df$gene_symbol,colnames(heatmap_input)),]
    
    #colnames(heatmap_input)
    #
    #df$gs_subcat<-str_to_title(tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",df$gs_subcat))))
    
    df$gs_subcat<-paste0(str_to_title(sub("(.*?)\\s(.*)","\\1",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",df$gs_subcat))))),
           " ",sub("(.*?)\\s(.*)","\\2",tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",df$gs_subcat)))))
    #Split into several spaces
    df<-df %>% mutate(gs_subcat= sapply(gs_subcat, function(x) paste(strwrap(x, 30), collapse = "\n")))
    
     ## creating dataframe for format selection
    #ask_annotation=data.frame()
    #color_list=c("YlOrBr","Greens","Blues","Oranges","YlGnBu","OrRd","PuBu","Reds")
    #color_list=rep(rownames(brewer.pal.info)[c(1,3,5,6,7,8,9,10,12,14,16,18,20,22,24,26)],times=2)
    #color_list=rep(c("Set1","Set3"),times=5)
    #class_list=str_to_title(tolower(str_sub(sub(".*?\\_","",c(unique(BBDD_temp[,sel_columnToFilterGroups])))))
    #class_list=c(unique(BBDD_temp[,sel_columnToFilterGroups]))
    
    #   ### asking defaults
    #   ask_annotation=rbind(ask_annotation,
    #                        cbind("Annotation parameter"=c(paste0(" annotation choose base color:"),
    #                                                       paste0(" annotation name:"),
    #                                                       paste0(" annotation color levels 1-9:")),
    #                              Choices=c(color_list[(j)],
    #                                       str_sub(class_list[(j)],1,50), 
    #                                        #str_sub(gsub("HALLMARK_|WP|REACTOME|KEEG","", class_list[(j)]),1,50),
    #                                       #str_to_title(tolower(str_sub(sub(".*?\\_","", class_list[(j)]),1,50))),
    #                                        "5,7,9"))
    #                        
    #   )
    #   
    #   
    # 
    # #print(lapply(df,unique)[-1])
    # 
    # ## select color for each option
    # fix(ask_annotation)
    # ask_annotation$Choices=as.character(ask_annotation$Choices)
    # 
    # 
    # ## creating color pallete and anotation for second clasification 
    # 
    # escala_x1=c(brewer.pal(9,ask_annotation[1,2])[as.numeric(str_split(ask_annotation[3,2],",",simplify = T)[1:length(unique(df[,1]))])])
    # names(escala_x1)=c(unique(df[,1])[order(unique(df[,1]))])
    # escala_x1[which(names(escala_x1)!=class_list[(j)])]="white"
    #  
    # Easier alternative
    aux_scale<-unlist(lapply(seq_len(length(ask_annotation)),function(i){
      c(brewer.pal(brewer.pal.info[ask_annotation[i],]$maxcolors,ask_annotation[i]))}))
    escala_x1<-c("white",aux_scale[j])
    names(escala_x1)=c(unique(df[,1])[order(unique(df[,1]))])
    #names(escala_x1)[2]<-str_to_title(tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",names(escala_x1)[2]))))
    
    ## saving annotation    
    anolist[[j]]<-rowAnnotation(Category = df[,1],
                          col = list(Category=escala_x1),
                          show_annotation_name = T, #False if not rows 
                          na_col = "white",
                          #name = ask_annotation[2,2],
                          name = names(escala_x1)[2],
                          #name = j,
                          #width=unit(5,"mm"),
                          # annotation_label=class_list[(j)]
                          annotation_label=as.character(j),
                          annotation_name_rot = 0,
                          annotation_name_side = "top"
                          )
    #anolist[[j]]@name<-str_to_title(tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",anolist[[j]]@name))))
    
    
    
    
  }

 
  
construct_annotation<-unlist(lapply(seq_len(length(anolist)), function(i){
    paste0("anolist[[",i,"]]")
  }))


  
# tt<-c("WP_CHOLESTEROL_BIOSYNTHESIS_DJKJKKDJFK_DFKKDJFKJF_fjkdjfk_dkjfkjdf")
# # tt<-c("Aniso (ddd)")
# aa<-str_to_title(tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",tt))))
#   #tolower(gsub("\\_"," ",sub("(.*?)\\_(.*)","\\2 (\\1)",tt)))
# #gsub(anolist[[1]]@name)

# paste(strwrap(aa, width=30),collapse = "\n")
# 
# aa<-str_split(aa, " ", n=3)
# paste0(aa[[1]][1]," ",aa[[1]][2],"\n",aa[[1]][3])
# 
# if(length(BBDD_sel_tto)==0){
#   ### FILTRAR LOS GRUPOS
#   sel_tto=selectOption(as.character(unique(gsub("(.*)\\D(\\d+$)","\\1",rownames(heatmap_input)))),name="toRemove",toRemove = T)
#   if(sum(!is.na(sel_tto))>0){
#     heatmap_input<-heatmap_input[gsub("(.*)\\D(\\d+$)","\\1",rownames(heatmap_input)) %in% sel_tto,]
#   }
#   rownames(heatmap_input)
# }

# if(length(auxVars)>0){
#   
#   ## checking ids are correct
#     gsub("(.*)\\D(\\d+$)","\\2",colnames(t(heatmap_input)))==
#     paste0(BBDD_temp$SAMPLE.ID,BBDD_temp$Batch.Number)
#   
#   heatmap_input=cbind(heatmap_input,(BBDD_temp[,auxVars]))
# }



#### HEATMAP SIMPLE #####
selection=selectOption(unique(gsub("(.*)\\D(\\d+$)","\\1",rownames(heatmap_input))),"toRemove",toRemove = T)
  
heatmap_input = dfc.dots[,-1]


heatmap_input = t(scale((heatmap_input)))

#BBDD_temp[BBDD_temp$gs_name %in% catlist,"gene_symbol"]

#Adjust number of rows to selected gene_name categories TODO: See what happened with other criteria  
heatmap_input<-heatmap_input[rownames(heatmap_input) %in% BBDD_temp[BBDD_temp$gs_name %in% catlist,"gene_symbol"],]

#heatmap_input=log2(t(normalizeTo1(dfc.dots[,-1][grep(paste0(paste0("^",selection),collapse="|"),rownames(dfc.dots[,-1])),])))
#heatmap_input = data.frame(t(heatmap_input))

##### promedio de log
# heatmap_input_bkp <- heatmap_input
# heatmap_input <- heatmap_input_bkp
# 
# heatmap_input <- as.data.frame(heatmap_input)
# heatmap_input$classes <- (gsub("_\\d+$","",rownames(heatmap_input)))
# 
# heatmap_input <- heatmap_input %>% group_by(classes) %>% summarise_if(is.numeric,mean) %>% as.data.frame()
# rownames(heatmap_input) <- heatmap_input$classes
# heatmap_input <- heatmap_input[,-1]

#### ordenar las clases

asklevels=c("Cambiar order clases?"="y")
fix(asklevels)

if(tolower(asklevels)=="y" | !exists("heatmap_levels")){
  levels=unique(gsub("(.*)\\D(\\d+$)","\\1",rownames(heatmap_input)))
  if(!exists("heatmap_levels")){
    heatmap_levels=as.data.frame(cbind(levels=levels,order=c(1:length(levels))),stringAsFactor=F)
  } else if(sum(heatmap_levels$levels %in% levels)!=length(levels)){
    heatmap_levels=as.data.frame(cbind(levels=levels,order=c(1:length(levels))),stringAsFactor=F)
  }
  fix(heatmap_levels)
}

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

#colnames(heatmap_input)=gsub("F4","F3",colnames(heatmap_input))
#cluster_within_group(heatmap_input)

# htmp_simple=draw(
#   Heatmap(heatmap_input,
#           name = "log2FC",
#           width=ncol(heatmap_input)*unit(15,"mm"),
#           height=nrow(heatmap_input)*unit(3,"mm"),
#           border_gp = gpar(col = "darkgrey", lty = 1,alpha=1),
#           rect_gp = gpar(col = "white", lwd = .2,lty = 1),
#           #column_title = "Heatmap rat?n", ## quita titulos del split
#           col = col_fun,
#           #ordenamiento
#           
#           #cluster_columns = cluster_within_group(heatmap_input, group),
#           #cluster_column_slices = F,
#           #cluster_row_slices = F,## tiene que ser false para poder ordenar los grupos 
#           #clustering_method_rows = "median",#centroid,median,mcquitty,average,complete,single,ward.D2,ward.D
#           #cluster_rows = function(m) as.dendrogram(diana(m)),
#           #cluster_columns = function(m) as.dendrogram(diana(m)),
#           #cluster_rows = T,
#           #clustering_distance_rows = "pearson",
#           #name = "mat", 
#           cluster_rows = diana,
#           cluster_columns=F,
#           clustering_method_columns="median",
#           #column_km = 3,
#           #nombres filas/columnas
#           show_row_dend  = F,
#           column_dend_reorder = T, #controla ordenamiento del dendrograma
#           show_row_names = T,
#           show_column_names = F, 
#           show_column_dend = F,
#           # parametros split
#           # row_split = factor(rownames(heatmap_input),
#           #                    levels=as.character(heatmap_levels$levels)[match(rownames(heatmap_levels),heatmap_levels$order)],ordered = T),
#           
#           #column_split = factor(gsub("(.*)\\D(\\d+$)","\\1",colnames(heatmap_input))),
#           # column_split = factor(gsub("(.*)\\D(\\d+$)","\\1",colnames(heatmap_input)),
#           #                       levels=as.character(heatmap_levels_group$levels)[match(rownames(heatmap_levels_group),heatmap_levels_group$order)],ordered = T),
#           column_split=factor(gsub("(.*)_(\\d+$)","\\1",colnames(heatmap_input)),
#                               levels=as.character(heatmap_levels_group$levels)[match(rownames(heatmap_levels_group),heatmap_levels_group$order)],ordered = T),
#           column_title_side = "top",
#           #column_title_gp = gpar(fontsize = 10), #Cajitas para los titulos
#           row_gap = unit(c(1), "mm"),
#           column_gap = unit(c(2), "mm"),
#           ### etiquetas de casos en filas
#           #right_annotation  = c(anolist[[1]],anolist[[2]],anolist[[3]],anolist[[4]],anolist[[5]],anolist[[6]],anolist[[7]],anolist[[8]]),
#           right_annotation  = eval(parse(text=paste0("c(",paste(construct_annotation,collapse=","),")"))), # to build annotation
#           column_names_rot = c(45),
#           column_title_rot = 0,
#           column_title_gp = gpar(fontsize = 15),
#           row_title_gp = gpar(fontsize = 15),
#           column_names_gp = gpar(fontsize = 10),
#           row_names_gp = gpar(fontsize = 10),
#           row_title_rot = 0
#           
#   )
# )



htmp<-Heatmap(heatmap_input,
          name = "log2FC",
          width=ncol(heatmap_input)*unit(15,"mm"),
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
          #right_annotation  = c(anolist[[1]],anolist[[2]],anolist[[3]],anolist[[4]],anolist[[5]],anolist[[6]],anolist[[7]],anolist[[8]]),
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

#draw(htmp, heatmap_legend_side="right", annotation_legend_side="top",
 #    legend_grouping = "original")

png(paste0(directory,"/Figures/HeatMap_aux_test_",as.numeric(Sys.time()),".png"),
    width = 60, height = 80, units = "cm", res = 200)
print(draw(htmp, heatmap_legend_side="bottom", annotation_legend_side="right",
           legend_grouping = "original")
)
dev.off()
browseURL(paste0(directory))







