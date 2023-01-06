
directory=get("directory",envir=.GlobalEnv)
print("Estas son las carpetas que hay ahora mismo creadas para analisis:")
print(matrix(gsub("(.*)\\/(.*)$","\\2",list.dirs(directory))))
carpeta_para_figuras=c("figurasQC")
names(carpeta_para_figuras) <-c("Introduzca un nombre para la carpeta donde se van a guardar las figuras generadas: ")
fix(carpeta_para_figuras)
#creamos directorio nuevo
dir.create(sprintf(paste0(directory,"/",carpeta_para_figuras)),showWarnings=F);
sub_folder<-paste0(directory,"/",carpeta_para_figuras)
print(paste("Directorio",sub_folder,"creado"))


#eindex = rsALL.samples[,c("index","NAME")]
#inputData=t_rsALL.concTOistd.i

##cargamos el archivo de QC
QC_file <-selectOption(list.files(directory)[grepl(".*(.csv)$",list.files(directory),ignore.case=TRUE)],"Select_QC_files")

QC_file<-file.path(directory,QC_file)
#cargamos
inputDataQC<-(read.csv(
  QC_file,
  header = T,
  sep = ",",
  dec = ".",
  stringsAsFactors = F,
  check.names = F,
  blank.lines.skip	= T
))



### seleccionamos la muestra a comparar como QC
select_qc=selectOption(unique(gsub("(\\D+)(\\d+$)",
                                          "\\1",
                                          gsub(
                                            "\\..*", "", rownames(inputData)
                                          ))),"Select_QC_SAMPLE")

select_qc_historical=selectOption(unique(gsub("\\#.*",
                                   "",inputDataQC$type)),"Select_QC_SAMPLE_historical")
 
inputData = inputData[grepl(select_qc, rownames(inputData)),]

## adaptamos el formato para mantener coherencia con la app
saveColOrderQC=colnames(inputData)
QC_method=c("ISTD/Cal")
names(QC_method)=c("Que analisis estas haciend? ISTD o Cal")
fix(QC_method)

inputData$type="SAMPLE"
inputData$Batch.number=rsSample[which(rsSample$Muestra %in% rownames(inputData)),"Tanda"]
inputData=inputData[,c("type","Batch.number",saveColOrderQC)]



##fix for QC historical file loaded

  
  if (sum(c("type", "Batch.number") %in% str_trim(colnames(inputDataQC)[1:2])) ==
      2) {
    ## comprobando que las especies coincidan 
    colnames(inputDataQC)[1:2] = c("type", "Batch.number")
    
    inputDataQC$type = gsub("\\#\\d+$", "", inputDataQC$type)
    eindex$NAME=gsub("_"," ",eindex$NAME)
    
    #inputDataQC=colnames(inputDataQC)[-(1:2)]
    #inputData=str_trim(eindex[grepl(paste0("^",str_trim(colnames(inputDataQC)[-(1:2)]),".*",collapse = "|") , str_trim(eindex$NAME)), "NAME"])
    ## convertimos en long format para poder analizar mas facilmente
    inputDataQC=reshape2::melt(inputDataQC, id.vars=c("type","Batch.number"))
    
    inputData=cbind(Sample.Name=rownames(inputData),inputData)
    colnames(inputData)[-c(1:3)]=str_trim(eindex[match(str_trim(colnames(inputData)[-(1:3)]) , str_trim(eindex$index)), "NAME"])
    
    inputData=reshape2::melt(inputData, id.vars=c("Sample.Name","type","Batch.number"))
    
   ## buscamos las especies del QC que pueden tener sufijos
    SpeciesMatch = inputDataQC %>% group_by(variable) %>% summarise(AnalysisMatch =
                                                                    paste0(unique(inputData[grepl(paste0("^", unique(variable), ".*"), (inputData$variable)), "variable"]), collapse = ","))
    ### por problema de columna factor
    SpeciesMatch$variable=as.character(SpeciesMatch$variable)
    
    ## corregimos los nombres de especies de analisis segun los nombres de especies en QC
    inputData = inputData %>% group_by(variable) %>% summarise(
      Sample.Name=Sample.Name,
      type = type,
      Batch.number = Batch.number,
      value = value,
      variable = as.character(unique(SpeciesMatch[grepl(paste0("^", unique(variable), ".*|\\,", unique(variable), ".*"),
                                                             (SpeciesMatch$AnalysisMatch)), "variable"]))
    )
    ## quitamos especies que no coinciden
    inputData=inputData[(inputData$variable)!="character(0)",]
    ## sumamos las especies duplicadas
    inputData=inputData%>%group_by(Sample.Name,variable,Batch.number)%>%summarise(variable=variable[1],type=type[1],value=sum(value,na.rm=T))
    
    
    
    inputDataQC = inputDataQC[as.character(inputDataQC$variable) %in% unique(inputData$variable),]
    
    inputDataQC = inputDataQC[grepl(paste0(select_qc_historical,"#",QC_method), inputDataQC$type), ]
    
  } 
  


### getting mean of actual experiment

dfc.analysis = inputData[,c("type", "Batch.number", "variable", "value")]
dfc.analysis = as.data.frame(
  dfc.analysis %>% group_by(variable, Batch.number) %>% summarise(
    type = type[1],
    Batch.number = Batch.number[1],
    variable = variable[1],
    value = my.median(value),
    sd = my.median(value)
  ),
  stringAsFactor = F
)
dfc.QC = as.data.frame(
  inputDataQC %>% group_by(variable) %>% summarise(sdu =
                                                     my.sd(value),value = my.median(value)),
  stringsAsFactors = F
)
dfc.QC$variable=as.character(dfc.QC$variable)
dfc.QC$value=as.numeric(dfc.QC$value)
### adaptando el archivo cargado

  
  ### multiplicating class mean for every species
  dfc.analysis_expanded = dfc.analysis
  dfc.analysis_expanded$consensus = (dfc.QC[match(dfc.analysis_expanded$variable,
                                                        dfc.QC$variable), "value"])
  
  
  ### getting SD
  dfc.analysis_expanded$sdu = (dfc.QC[match(dfc.analysis_expanded$variable,
                                            dfc.QC$variable), "sdu"])
  
  #colnames(dfc.srm_fus_sdu) = c("type", "Batch.number", "variable", "sdu")
  

## visualization
dfc.analysis_expanded$Zscore <-
  (dfc.analysis_expanded$consensus - dfc.analysis_expanded$value) / dfc.analysis_expanded$sdu
dfc.analysis_expanded$Zscore = ifelse(dfc.analysis_expanded$Zscore > 10,
                            10,
                            ifelse(dfc.analysis_expanded$Zscore < (-10), -10, dfc.analysis_expanded$Zscore))
dfc.analysis_expanded$SDnorm <- dfc.analysis_expanded$sd / dfc.analysis_expanded$sdu


dfc.analysis_expanded$outlier_exactitud <-
    ifelse(abs(dfc.analysis_expanded$Zscore) > 3, 0.9, ifelse(abs(dfc.analysis_expanded$SDnorm) > 2, 0.9, 0.5))


### seleccionar clases a mostrar
  #hideClass = selectOption(unique( gsub("_.*","",dfc.analysis_expanded$variable)),"Select_class_to_show")
  #dfc.analysis_expanded = dfc.analysis_expanded[grepl(hideClass, dfc.analysis_expanded$variable), ]



dfc.analysis_expanded$Class = gsub("(^\\S+)\\s.*", "\\1", dfc.analysis_expanded$variable)

dfc.analysis_expanded$Batch.number=as.character(dfc.analysis_expanded$Batch.number)
#pintamos las graficas
 
  dfc.analysis_expanded$variable=gsub("(^\\w*)(.*)", "\\2",dfc.analysis_expanded$variable)
  dfc.analysis_expanded[is.na(dfc.analysis_expanded$variable),"variable"]=dfc.analysis_expanded[is.na(dfc.analysis_expanded$variable),"Class"]
  
  palete = colorRampPalette(c("red","brown"))(length(unique(dfc.analysis_expanded$Batch.number)))
  pd <- position_dodge(1)
  
  
  for(cls in dfc.analysis_expanded$Class){
    sub_dfc.srm_fus=dfc.analysis_expanded[which(dfc.analysis_expanded$Class %in% cls),]
    
    myPng(paste(paste(cls,collapse="_"),"_QC_monitoring",sep=""),  16,9, 200,dir=sub_folder);
    
  
  p <-
    ggplot(sub_dfc.srm_fus, aes(x =variable, y = Zscore)) +

    geom_rect(
      xmin = 0,
      xmax = Inf,
      ymin = -2,
      ymax = 2,
      fill = "#c0ffc0"
    ) +
    geom_rect(
      xmin = 0,
      xmax = Inf,
      ymin = 2,
      ymax = 3,
      fill = "#ffc0c0"
    ) +
    geom_rect(
      xmin = 0,
      xmax = Inf,
      ymin = -2,
      ymax = -3,
      fill = "#ffc0c0"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "black",
      size = 0.5
    ) +
    geom_hline(
      yintercept = 2,
      linetype = "dashed",
      colour = "green",
      size = 0.5
    ) +
    geom_hline(
      yintercept = -2,
      linetype = "dashed",
      colour = "green",
      size = 0.5
    ) +
    geom_hline(
      yintercept = 3,
      linetype = "dashed",
      colour = "red",
      size = 0.5
    ) +
    geom_hline(
      yintercept = -3,
      linetype = "dashed",
      colour = "red",
      size = 0.5
    ) +
    geom_text(
      x = Inf,
      y = 2.1,
      label = "2-sd",
      size = 3,
      color = "black",
      hjust = 1.2
    ) +
    geom_text(
      x = Inf,
      y = 3.1,
      label = "3-sd",
      size = 3,
      color = "black",
      hjust = 1.2
    ) +

    
    guides(alpha = "none") +
    geom_point(
      aes(size = Batch.number,
          shape = Batch.number),
      colour = "black",
      size = 2,
      stroke = 1,
      alpha = .5,
      show.legend = F
    ) +
    geom_point(
      size = 2,
      aes(
        color = Batch.number,
        fill = Batch.number,
        size = Batch.number,
        shape = Batch.number
        # alpha = outlier_exactitud
      )) +
    scale_shape_manual(values = c(22:26)) +
    scale_fill_manual(values =  palete[1:length(unique(sub_dfc.srm_fus$Batch.number))]) +
    scale_color_manual(values =  palete[1:length(unique(sub_dfc.srm_fus$Batch.number))])

  
  
  # 
  # if (input$ISTD_QC_Plot_sdu) {
  #   p = p + geom_errorbar(
  #     aes(
  #       ymin = Zscore - SDnorm,
  #       ymax = Zscore + SDnorm,
  #       alpha = outlier_exactitud,
  #       color = Batch.number
  #     ),
  #     #alpha = 0.5,
  #     width = 1,
  #     show.legend = F,
  #     size = 1
  #   )
  # }
  #geom_text(x=Inf,y=-3.1,label="3-sd",size=3,color="grey",hjust=1.2)
  
  
  #myPng("Pecision_Coverage_species_inLOQ", 9.57,6, 300);
  
  p = p + facet_wrap(~ Class,
                     scale = "free",
                     ncol = 4) +
    labs(
      title = paste0(
        "Between batch monitoring (",
        select_qc,
        ")"
      )
    ) +
    xlab("") +
    ylab("Coverage Equivalent (k-eq)") + #ojo con unidades
    #ylab(expression(paste("% of total")))+ #ojo con unidades
    theme_bw(15) +
    #sISTDe_y_continuous(trans='pseudo_log')+
    theme(
      #legend.position=c(.63, .9),
      legend.position = "top",
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      # switch off major gridlines
      panel.grid.minor = element_blank(),
      # switch off minor gridlines
      panel.border = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.line = element_line(size = 0.5),
      axis.ticks =  element_line(size = 0.3),
      axis.text.y = element_text(
        size = 14,
        angle = 0,
        face = "plain"
      ),
      axis.text.x = element_text(
        size = 10,
        angle = 90,
        hjust = 1.0,
        vjust = 0.5,
        face = "plain"
      ),
      #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
      #panel.border = element_border(c("left","bottom")),
      legend.title = element_blank(),
      # switch off the legend title
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.0, "lines"),
      axis.title.y = element_text(size = 16)
    )
  print(p)
  
  pushViewport(viewport(y=0.95,x=0.85,height=.05))
  title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
  grid.table(title)
  
  dev.off()
  }

