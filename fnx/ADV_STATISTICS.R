
summaryRES_PACK<-function(dfc_obj, digitos=2,varTOcast,estat="media_se"){
  #funcion recibe un objeto dfc
  #digitos=2 numero de digitos de redondeo
  #varTOcast nombre del la variable que agrupa los casos
  #estat, por defecto la media_sd, mediana_ci estat="mediana", media_se, estat=media_se
  #dfc_obj<-dfc.fibro;
  #digitos<-2;
  dfc_obj$msd<-paste0(paste(round(dfc_obj$value,digits=digitos),round(dfc_obj$sd,digits=digitos),sep=" \u00b1 ")," (",dfc_obj$N,")")
  dfc_obj$medci<-paste0(paste(round(dfc_obj$median,digits=digitos),round(dfc_obj$ci,digits=digitos), sep=" \u00b1 ")," (",dfc_obj$N,")")
  dfc_obj$mse<-paste0(paste(round(dfc_obj$value,digits=digitos),round(dfc_obj$se,digits=digitos), sep=" \u00b1 ")," (",dfc_obj$N,")")
  
  #dfc_obj$NAME<-rsALL.i$NAME[match(dfc.fc$variable,rsALL.i$id_spec)]
  #dfc.fc$id_spec<-rsALL.i$id_spec[match(dfc.fc$variable,rsALL.i$id_spec)]
  rsdfc_obj<-reshape2::dcast(dfc_obj, paste("variable", "~",varTOcast),value.var="msd")
  if(estat=="mediana"){
    rsdfc_obj<-reshape2::dcast(dfc_obj, paste("variable", "~",varTOcast),value.var="medci") }
  else if (estat=="media_se"){
    rsdfc_obj<-reshape2::dcast(dfc_obj, paste("variable", "~",varTOcast),value.var="mse") }
  else {
    rsdfc_obj<-reshape2::dcast(dfc_obj, paste("variable", "~",varTOcast),value.var="mse")
  }
  return(rsdfc_obj) 
}


summaryRES_ttest<-function(t_aux,varTOgroup,varTOind,varTOtest=NULL,pareado=F,alfa=0.05){
  
  #  funcion recibe un objeto t_rsALL con columnas casos-individuos (varTOind) y grupos (varTOgroup)
  #  grupos varTOgroup 
  #  individuos (en caso que se quiera test pareados) varTOind
  #  variables para realizar estadistica varTOtest
  #  falta inclusuion de test pareados
  
  #t_aux<-t_rsRESUM.esp
  #varTOgroup<-c("F_fibrosis")
  #varTOind<-c("fusion")
  #varTOtest<-measure_vars[1:3]
  
  t_aux<-t_aux[,c(varTOind,varTOgroup,varTOtest)];
  t_aux<-t_aux[order(t_aux[[varTOgroup]],t_aux[[varTOind]]),] 
  if(pareado==T){
    print(varTOtest)}  ##TODO test pareados
  t_aux<-t_aux[,-c(1)] #elimina la condicion de individuo
  t_aux[,c(1)]<-as.character(t_aux[,c(1)])
  colnames(t_aux)[1]<-"tto"
  
  ##sustituir nombres de variables por nombres provisionales para evitar fallos en llamada a superttest
  string<-list();
  for (i in 1:length(letters)){string[[i]]<-paste(LETTERS,letters[i],sep="")}
  string<-unlist(string)
  
  colnames(t_aux)<-string[1:length(t_aux)]
  colnames(t_aux)[1]<-"tto"
  #unique(t_aux$tto)
  
  ###cargar las combinaciones posibles de comparaciones###
  dfx<-expand.grid(unique(t_aux$tto),unique(t_aux$tto),stringsAsFactors = F)
  dfx<-dfx[dfx$Var1 <= dfx$Var2,] 
  dfx<-dfx[!dfx$Var1==dfx$Var2,]
  dfx<-dfx[order(dfx[,1]),]
  t_test<-list()   
  
  for(i in 1:length(dfx[,1])){
    t_test[[i]]<-super_T_test(dat_test=t_aux,vector1=dfx[i,1],vector2=dfx[i,2],rep=F)
  }
  output <- t(matrix(unlist(t_test),ncol=length(dfx[,1]),nrow=length(varTOtest)))
  colnames(output)<-varTOtest
  output<-cbind(dfx,output)
  output[,2]<-paste(output[,1],output[,2],sep=" to ")
  output<-output[,-c(1)]
  nam<-output$Var2
  output<-t(output[,-c(1)])
  colnames(output)<-nam
  output[which(output>=alfa)]<-NA; 
  output<-round(output,digits=3);
  #contabiliza n comparaciones significativas
  if (dim(output)[2]>1){
    N_sig<-apply(output[,which(!colnames(output)=="variable")],1,function(x){length(x[!is.na(x)])}); #exluding QC value, #listado datos que son 0 o NA
  } else{
    N_sig<-as.numeric(!is.na(output[,which(!colnames(output)=="variable")]))
  }
  output<-cbind(output,N_sig)
  return(output)     
}


multTTest<-function(df,set1,set2,group.var,repeated){
  #if("dplyr" %in% (.packages())){
  # detach("package:reshape", unload=TRUE)
  #detach("package:reshape2", unload=TRUE)
  #detach("package:dplyr", unload=TRUE)
  #detach("package:plyr", unload=TRUE)
  #} 
  #library(plyr)
  #library(reshape)
  #library(reshape2)
  #library(dplyr)
  
  ###Hay problemas de compatibilidad entre plyr y dplyr que no se resover
  formula<-as.formula(paste(group.var,"tto",sep="~"))
  #df$tto<-as.character(df$tto)
  dt_res <-
    data.frame(expand.grid(set1,set2)) %>%  # create combinations 
    dplyr::mutate(test_id = row_number()) %>%    # create a test id
    group_by(test_id) %>%              # group by test id, so everything from now on is performed for each test separately
    do({x_temp = df[(df$tto == .$Var1 | df$tto == .$Var2),]    # for each test id keep groups of interest
    x_temp = data.frame(x_temp)}) %>%
    #do(test = t.test(data~Group, data=.)) 
    do(test = t.test(formula, data=.,na.action=na.omit,paired = repeated)) 
  
  #detach("package:dplyr", unload=TRUE)
  
  return(dt_res) 
}


super_T_test<-function(dat_test,vector1,vector2,rep=F){
  
  tvector<-colnames(dat_test)[-1]
  tt_plus<-lapply(tvector,
                  FUN=function(x) multTTest(
                    df=dat_test,
                    set1=vector1,
                    set2=vector2,
                    group.var=x,
                    repeated=rep
                  )
  )
  
  p.value<-lapply(1:length(tt_plus),function(x){
    sapply(tt_plus[[x]]$test, getElement, name = "p.value")
  })
  p.names_comp<-lapply(1:length(tt_plus),function(x){
    #aux2<-list();
    lapply(1:length(tt_plus[[x]]$test_id),function(y){
      aux<-gsub("mean in group ","",names(tt_plus[[x]]$test[[y]]$estimate));
      aux<-paste(aux[1],aux[2],sep=" to ");
      return(aux)
    })
  })
  p.values_matrix<-matrix(unlist(p.value), ncol =length(tt_plus), byrow = F)
  colnames(p.values_matrix)<-tvector
  rownames(p.values_matrix)<-unlist(p.names_comp[[1]])
  p.values_matrix<-p.values_matrix[order(rownames(p.values_matrix)),]
  
  return(p.values_matrix)
}

summaryRESSTAT <-
  function(t_frame,
           case_vars,
           measure_vars,
           varTOcast,
           digitos,
           estat,
           sum = F,
           units="value") {
    #recibe un objeto t_frame (data.frame)
    #funcion recibe un objeto dfc
    #digitos=2 numero de digitos de redondeo
    #varTOcast nombre del la variable que agrupa los casos
    #estat, por defecto la media_sd, mediana_ci estat="mediana", media_se, estat=media_se
    #dfc_obj<-dfc.fibro;
    #digitos<-2;
    ### devuelve una lista de objetos (t_frame,m_frame,dfc,table) para poder seguir trabajando
    
    # t_frame=inputDataCalraw
    # case_vars="GROUP"
    # measure_vars= colnames(inputDataCalraw)[-c(1, 2, (ncol(inputDataCalraw) - 1), ncol(inputDataCalraw))]
    #
    # varTOcast = "GROUP"
    # digitos = 3
    # estat = "media_se"
    # sum = T
    # save objects in current environment
    # save(list = ls(), file = "shiny_summaryRESSTAT_qc2.Rdata", envir = environment())
    # load("shiny_summaryRESSTAT_qc2.Rdata")

    # browser()
    # detach("package:plyr",unload=T)
    t_frame$ID=paste0(t_frame$GRUPO,"_",t_frame$ID)
    m_frame <-
      reshape2::melt(t_frame, id = case_vars, measure = measure_vars)
    m_frame[] <- lapply(m_frame, as.character)
    m_frame$value=as.numeric(m_frame$value)
    
    if (sum == T) {
      # browser()
      m_frame$class = gsub("\\s.*", "", m_frame$variable)
      
      m_frame = ddply(m_frame, .(ID, class), colwise(myfun2, c("value")))
      colnames(m_frame) = c("GRUPO", "variable", "value")
      t_frame = reshape2::dcast(m_frame, GRUPO ~ variable, value.var = "value")
      m_frame$GRUPO = gsub("(\\D+)(\\d+$)", "\\1", m_frame$GRUPO)
      
    } 
    
    
    #t_frame = t_frame[, -c((ncol(t_frame) - 2):ncol(t_frame))]
    dfc <-
      summarySE(
        m_frame,
        measurevar = "value",
        groupvars = c(varTOcast, "variable"),
        na.rm = T
      )
    
    
    
    table <-
      summaryRES_PACK(dfc,
                      varTOcast = varTOcast,
                      estat = estat,
                      digitos = digitos)
    my_groups <- unique(dfc[, 1])
    dfx <- expand.grid(my_groups, my_groups, stringsAsFactors = T)
    dfx <- dfx[as.character(dfx$Var1) <= as.character(dfx$Var2),]
    dfx <- dfx[!dfx$Var1 == dfx$Var2,]
    dfx <- dfx[order(dfx[, 1]),]
    my_comparisons <-
      mapply(c,
             as.character(dfx$Var1),
             as.character(dfx$Var2),
             SIMPLIFY = F)
    symnum_list <-
      list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        symbols = c("<0.0001", "<0.001", "<0.01", "<0.05", "ns")
      )
    ##ANOVA and t.test

    ANOVA_sig <-
      compare_means(
        as.formula(paste0('value', '~', varTOcast)),
        data = m_frame,
        method = "anova",
        group.by = c("variable"),
        symnum.args = symnum_list
      )
    check_replicates = table(t_frame$GRUPO)
    if (sum(check_replicates == 1) < length(check_replicates)) {
      check_replicates=check_replicates[!(check_replicates==1)]
      m_frame=m_frame[which(m_frame$GRUPO %in% names(check_replicates)),]
      if(length(unique(m_frame$value))<length(m_frame$value)){
        
        table <-
          join(table, ANOVA_sig[, c("variable", "p.signif")], by = "variable", type =
                 "left") ##fusionar
      } else {
        
        
        t_test_sig <-
          compare_means(
            as.formula(paste0('value', '~', varTOcast)),
            data = m_frame,
            comparisons = my_comparisons,
            method = "t.test",
            group.by = c("variable"),
            symnum.args = symnum_list
          ) #+ # Add pairwise comparisons p-value
        t_test_sig$compara <-
          paste(t_test_sig$group1, t_test_sig$group2, sep = "_vs_")
        t_test_sig_obj <-
          reshape2::dcast(t_test_sig, variable ~ compara, value.var = "p.signif")
        #contabiliza n comparaciones significativas
        if (dim(t_test_sig_obj)[2] > 2) {
          N_sig <-
            apply(t_test_sig_obj[, which(!colnames(t_test_sig_obj) == "variable")], 1, function(x) {
              length(x[!is.na(x) &
                         !x == "ns"])
            })
          #exluding QC value, #listado datos que son 0 o NA
        } else{
          N_sig <-
            as.numeric(!is.na(t_test_sig_obj[, which(!colnames(t_test_sig_obj) == "variable")]) &
                         !t_test_sig_obj[, which(!colnames(t_test_sig_obj) == "variable")] == "ns")
        }
        t_test_sig_obj <- cbind(t_test_sig_obj, N_sig)
        table <-
          join(table, ANOVA_sig[, c("variable", "p.signif")], by = "variable", type =
                 "left") ##fusionar
        table <-
          join(table, t_test_sig_obj, by = "variable", type = "left") ##fusionar
      }
    } else {
      table <-
        join(table, ANOVA_sig[, c("variable", "p.signif")], by = "variable", type =
               "left") ##fusionar
    }
    if(sum==F){
      
      t_frame=t_frame[,-1]
      m_frame <-
        reshape2::melt(cbind(Sample.Name=rownames(t_frame),t_frame), id = c("ID"), measure = measure_vars)
      m_frame[] <- lapply(m_frame, as.character)
      m_frame$value=as.numeric(m_frame$value)
      colnames(m_frame)=c("Sample.Name","Lipid.Name",units)
    } else {
      m_frame <-
        reshape2::melt(t_frame, id = "GRUPO", measure = colnames(t_frame)[-1])
      m_frame[] <- lapply(m_frame, as.character)
      m_frame$value=as.numeric(m_frame$value)
      colnames(m_frame)=c("Sample.Name","Lipid.Name",units)
    }
    colnames(dfc)[1:4]=c("GRUPO","Lipid.Name","N",units)
    colnames(table)[1]=c("Lipid.Name")
    if(sum==F){
      reorder_wide=match(c("GRUPO","ID","Batch.number"),colnames(t_frame))
      t_frame=t_frame[,c(colnames(t_frame)[reorder_wide],colnames(t_frame)[-c(reorder_wide)])]
    }
    
    obj <-
      list(
        "wideData" = t_frame,
        "longData" = m_frame,
        "SummaryStats" = dfc,
        "ANOVA" = table
      )
    
    return(obj)
  }

expandRows <- function(dataset, count, count.is.col = TRUE) {
  #expande filas en funcion de cuenta
  # expandRows(mydf, 2, count.is.col=FALSE)
  # expandRows(mydf, "frequency")
  # expandRows(mydf, c(1, 2, 1, 0, 2), count.is.col=FALSE)
  if (!isTRUE(count.is.col)) {
    if (length(count) == 1) {
      dataset[rep(rownames(dataset), each = count), ]
    } else {
      if (length(count) != nrow(dataset)) {
        stop("Expand vector does not match number of rows in data.frame")
      }
      dataset[rep(rownames(dataset), count), ]
    }
  } else {
    dataset[rep(rownames(dataset), dataset[[count]]), 
            setdiff(names(dataset), names(dataset[count]))]
  }
}

QC <- function(t_esp,
               t_esp.i,
               eindex,
               units_to_use="pmol/ug",
               fixed_scale = 0
) {
  ##recibe una tablas y calcula la calidad de la cuantificacion
  ## La tabla debe incluir SAMPLEid = (identificacion de la muestra)
  ##                       tanda = si hay varias tandas de analisis
  ##                       type  = si es control (QC) o muestra (SAMPLE)
  ## Las imputaciones se realizan seg?n rownames que tiene que tener la estructura (GRUPO_sAMPLEid)
  # t_esp<-t_rsALL.concTOcal
  # t_esp.i<-t_rsALL.concTOcal[,-c(1:2)]
  # eindex
  # units_to_use = "pmol/ug"
  # # fixed_scale=palete
  # browser()
  # save(list = ls(), file = "shiny_QC.Rdata", envir = environment())
  # load("shiny_QC.Rdata")
  library(ade4)
  t_esp.QC <- t_esp[which(t_esp$type %in% c("QC")),]
  selectedQC = unique(gsub("(\\D+)(\\d+$)", "\\1", gsub("\\..*", "", rownames(t_esp.QC))))
  t_esp.QC <- cbind(SAMPLEid = rownames(t_esp.QC), t_esp.QC)
  
  colnames(t_esp.QC)
  m_esp <-
    reshape2::melt(t_esp.QC, id = colnames(t_esp.QC)[1:3], na.rm = F)
  m_esp$variable = as.character(m_esp$variable)
  
  m_esp = m_esp[!is.na(m_esp$value), ]
  ### CV intra-dia
  dfc.qc.es <-
    summarySE(
      m_esp,
      measurevar = "value",
      groupvars = c("variable", "tanda"),
      notification = F,
      na.rm = T
    )
  
  dfc.qc.cvi <-
    summarySE(
      dfc.qc.es,
      measurevar = "cv",
      groupvars = c("variable"),
      na.rm = T
    )
  
  ### CV inter-dia
  dfc.qc.es.inter <-
    summarySE(
      dfc.qc.es,
      measurevar = "value",
      groupvars = c("variable"),
      na.rm = T
    )
  
  #### CV total
  dfc.qc.es.total <-
    summarySE(
      m_esp,
      measurevar = "value",
      groupvars = c("variable"),
      na.rm = T
    )
  
  ######
  dfc.qc.es.total <-
    cbind(dfc.qc.es.total,
          cv_inter = dfc.qc.es.inter$cv,
          cv_intra = dfc.qc.cvi$cv)
  ######
  dfc.qc.es.total <-
    cbind(eindex[match(dfc.qc.es.total$variable, eindex$index),], dfc.qc.es.total)
  
  ## tabla resumen
  dfc.qc.es.total = dfc.qc.es.total[match(colnames(t_esp.QC[, -c(1:3)]), dfc.qc.es.total$variable), ]
  dfc.qc.es.inter = dfc.qc.es.inter[match(colnames(t_esp.QC[, -c(1:3)]), dfc.qc.es.inter$variable), ]
  dfc.qc.es = dfc.qc.es[match(colnames(t_esp.QC[, -c(1:3)]), dfc.qc.es$variable), ]
  dfc.qc.cvi = dfc.qc.cvi[match(colnames(t_esp.QC[, -c(1:3)]), dfc.qc.cvi$variable), ]
  
  
  #dfc.qc.es.total$Region<-sub("cyl","",dfc.qc.es.total$Region)
  ###Marcar especies con CV>20 y CV<20
  #dfc.qc.es.total_exc<-dfc.qc.es.total[which(dfc.qc.es.total$cv>=20|dfc.qc.es.total$cv_intra>=20|dfc.qc.es.total$cv_inter>=20),]
  #
  m_QC.esp.total <-
    reshape2::melt(dfc.qc.es.total[, c("Region", "cv", "cv_inter", "cv_intra")], id =
                     c("Region"))
  
  
  ###### HISTOGRAM M#####
  # Histogram on a Continuous (Numeric) Variable
  theme_set(theme_classic())
  
  #escala = distinctColorPalette(k = length(unique(m_QC.esp.total$Region)))
  if (fixed_scale[1] == 0) {
    escala = distinctColorPalette(k = length(unique(m_QC.esp.total$Region)))
  } else {
    escala = fixed_scale
  }
  escala = escala[1:length(unique(m_QC.esp.total$Region))]
  #cvto_plot<-"cv"
  m_QC.esp.total$variable = revalue(
    m_QC.esp.total$variable,
    c(
      "cv" = "RSD total %",
      "cv_inter" = "RSD between %",
      "cv_intra" = "RSD within %"
    )
  )
  
  g <-
    ggplot(m_QC.esp.total, aes(x = value, fill = Region)) + scale_fill_manual(values =  escala)
  
  
  #myPng("png_Hist_CV_especies_to_IS_all", 9.57, 5, 300);
  graf_cv_especies <- g + geom_histogram(
    aes(fill = Region),
    binwidth = 5,
    #bins=5,
    col = "black",
    size = .1,
    boundary = 0
  ) +
    labs(
      title = paste0("Histogram species (", selectedQC, ")"),
      subtitle = paste0(length(rownames(dfc.qc.es.total)), " species in ", length(unique(
        dfc.qc.es.total$Region
      )), " Regiones"),
      x = "",
      y = "N"
    ) +
    theme_bw() +
    facet_wrap( ~ variable,
                scales = "free",
                nrow = 1,
                ncol = 3) +
    #scale_y_continuous(expand = c(0, 0), limits = c(0, 500)) +
    # scale_x_continuous(breaks=0:10*5,expand=c(0,0), limits = c(0, max(m_QC.esp.total$value)+1))+
    scale_x_continuous(
      breaks = 0:10 * 5,
      expand = c(0, 0),
      limits = c(0, 35)
    ) +
    
    theme(
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = 0.5,
        face = "plain",
        color=rep(c("black","white"),times=7)
      ),
      axis.text.y = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = .5,
        face = "plain"
      ),
      #plot.margin = margin(2, 2, 2, 2, "cm"),
      strip.text.x = element_text(size = 12),
      #legend.position=c(.10, .98),
      #legend.title = element_blank(),
      legend.position = "none",
      axis.title.y = element_text(
        size = rel(1.3),
        angle = 0,
        vjust = 0.5
      ),
      axis.title.x = element_text(size = rel(1.3), angle = 0),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      # switch off major gridlines
      panel.grid.minor = element_blank(),
      # switch off minor gridlines
      panel.border = element_blank(),
      axis.line = element_line(size = 0.5),
      #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
      #panel.border = element_border(c("left","bottom")),
      legend.title = element_blank(),
      # switch off the legend title
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.0, "lines")
      
    )
  
  #dev.off()
  
  ###RESIDUOS
  ## obtengo la diferencia 100*((y-<y>)/<y>)
  ## correction del orden de columnas
  
  
  
  
  t_esp.QC_resid <- t(t_esp.QC[, -c(1:3)]) - dfc.qc.es.total$value
  t_esp.QC_resid <- t(t_esp.QC_resid / dfc.qc.es.total$value) * 100
  #esp.QC_resid<-t(t_esp.QC_resid)
  m_esp.QC_resid <- reshape2::melt(t_esp.QC_resid)
  colnames(m_esp.QC_resid)[3] <- "resid"
  m_esp.QC_resid = m_esp.QC_resid[!is.na(m_esp.QC_resid$resid), ]
  if (sum(m_esp.QC_resid$Var1 != m_esp$SAMPLEid) == 0) {
    m_esp <- cbind(m_esp, resid = m_esp.QC_resid$resid)
  } else{
    stop((
      "NO PUEDO HACER LOS RESIDUOS REVISA QUE TIENES LAS MISMAS COLUMNAS!!!"
    ))
  }
  #
  # escala <-
  #   c(
  #     '#e41a1c',
  #     '#377eb8',
  #     '#4daf4a',
  #     '#984ea3',
  #     '#fee090',
  #     '#ffff33',
  #     '#a65628',
  #     '#f781bf'
  #   )
  # escala <- escala[c(1, 5, 2, 7, 3, 6, 4, 8)]
  #escala = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  #escala=c("black","blue","red","green")
  #escala[index(escala)+1]
  if (fixed_scale[1] == 0) {
    escala_batch = distinctColorPalette(k = length(unique(m_esp$tanda)))
  } else {
    escala_batch = fixed_scale
  }
  escala_batch = escala_batch[1:length(unique(m_esp$tanda))]
  
  resid_1 <- ggplot(m_esp, aes(x = value, y = resid)) +
    geom_hline(yintercept = 0,
               color = "grey",
               linetype = 'dashed') +
    geom_hline(yintercept = 20,
               color = "grey",
               linetype = 'dashed') +
    geom_hline(yintercept = -20,
               color = "grey",
               linetype = 'dashed') +
    #geom_vline(xintercept=0,color="grey",linetype = 'dashed')+
    geom_point(
      aes(
        shape = as.character(tanda),
        fill = as.character(tanda)
      ),
      size = 2,
      stroke = 0.5,
      alpha = .80
    ) +
    scale_shape_manual(values = c(22:26)) +
    scale_fill_manual(values = escala_batch) +  ##para tener los mismos ticks
    scale_y_continuous(breaks = c(-(1:5) * 20, (1:5) * 20), limits = c(-50, 50)) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x)
        10 ^ x),
      labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                   .x))
    ) +
    labs(
      title = paste0("Residual plot all species (", selectedQC, ")"),
      x = units_to_use,
      y =  expression(paste("100 * (x - ", mu, " ) / ", mu))
    ) +
    theme(
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = 0.5,
        face = "plain"
      ),
      axis.text.y = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = .5,
        face = "plain"
      ),
      strip.text.x = element_text(size = 12),
      legend.position = "top",
      axis.title.y = element_text(
        size = rel(1.3),
        angle = 90,
        vjust = 0.0
      ),
      axis.title.x = element_text(size = rel(1.3), angle = 0),
      strip.background = element_blank(),
      axis.line = element_line(size = 0.5),
      legend.title = element_blank(),
      # switch off the legend title
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.0, "lines"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      )
    ) + annotation_logticks(sides = "tb")
  #resid_1
  ### Efecto en SUMAS ####
  
  #chequeamos si se han eliminado casos
  esp.i <- data.frame(t(t_esp.i))
  #esp.i<-cbind(Region=gsub("(.*)\\_(.*)","\\1",rownames(esp.i)),esp.i)
  esp.i <-
    cbind(Region = eindex[match(rownames(esp.i), eindex$index), "Region"], esp.i)
  
  
  SUM.i <-
    ddply(esp.i, .(Region), colwise(myfun2, colnames(esp.i)[-c(1)]))
  #Store each group sum
  t_SUM.i <- data.frame(t(SUM.i[, -c(1)]))
  #traspose
  colnames(t_SUM.i) <- SUM.i$Region
  
  #chequeamos si se han eliminado casos
  
  if (identical(rownames(t_esp), rownames(t_esp.i))) {
    #print("Las filas son identicas no se eliminado muestras!!!")
    t_SUM.i <-
      cbind(t_esp[rownames(t_esp.i), c("type", "tanda")], t_SUM.i)
  } else{
    #print(paste0("se elimino la muestra ", rownames(t_esp)[which(!rownames(t_esp) %in%
    #                                                               rownames(t_esp.i))]))
    t_SUM.i <-
      cbind(t_esp[rownames(t_esp.i), c("type", "tanda")], t_SUM.i)
  }
  
  #CV sumas
  t_SUM.QC <- t_SUM.i[which(t_SUM.i$type %in% c("QC")),]
  m_SUM.QC <-
    reshape2::melt(t_SUM.QC, id = c("tanda", "type"))
  m_SUM.QC = m_SUM.QC[(m_SUM.QC$value != 0), ]
  ### CV intra-dia
  dfc.qc.sum <-
    summarySE(
      m_SUM.QC,
      measurevar = "value",
      groupvars = c("variable", "tanda"),
      na.rm = T
    )
  
  dfc.qc.sum.cvi <-
    summarySE(
      dfc.qc.sum,
      measurevar = "cv",
      groupvars = c("variable"),
      na.rm = T
    )
  
  ### CV inter-dia
  dfc.qc.sum.inter <-
    summarySE(
      dfc.qc.sum,
      measurevar = "value",
      groupvars = c("variable"),
      na.rm = T
    )
  
  #### CV total
  dfc.qc.sum.total <-
    summarySE(
      m_SUM.QC,
      measurevar = "value",
      groupvars = c("variable"),
      na.rm = T
    )
  
  ######
  dfc.qc.sum.total <-
    cbind(dfc.qc.sum.total,
          cv_inter = dfc.qc.sum.inter$cv,
          cv_intra = dfc.qc.sum.cvi$cv)
  
  clase.set <- dfc.qc.sum.total[, c(1, 7, 10, 11)]
  orden <- order(clase.set$cv, decreasing = T)
  clase.set$variable <-
    factor(clase.set$variable,
           levels = unique(clase.set$variable)[orden],
           ordered = T)
  clase.set <- clase.set[order(clase.set$variable),]
  clase.set$paleta <- escala
  #rep(c(brewer.pal(12, "Paired"), "white"), times = 3)[1:length(clase.set$variable)]
  clase.set$paleta <-
    factor(clase.set$paleta,
           levels = unique(clase.set$paleta)[orden],
           ordered = T)
  
  ##AGRUPADO###
  m_clase.set <-
    reshape2::melt(clase.set, id.vars = c("variable", "paleta"))
  colnames(m_clase.set)[1] <- "clase"
  #myPng("png_Hist_CV_sums_to_IS_inLTR", 9.57, 5, 300);
  graf_cv_clases <-
    ggplot(m_clase.set, aes(variable, value, fill = clase)) +
    geom_bar(
      stat = "identity",
      position = "dodge",
      width = .7,
      size = 0.1,
      col = "black"
    ) +
    #  geom_text(aes(label=clase), position=position_dodge(width=0.9),angle=45, vjust=-0.25)+
    #geom_bar(aes(y = cv_intra, fill = paleta)) +
    scale_fill_manual(values = as.character(levels(m_clase.set$paleta)),
                      breaks = as.character(m_clase.set$clase)) +
    labs(
      title = paste0("RSD class sum (", selectedQC, ")"),
      #subtitle="total (cv), within-day (cv_intra), between-day (cv_inter)",
      x = "",
      y = "%"
    ) +
    scale_y_continuous(
      breaks = 0:5 * 10,
      expand = c(0, 0),
      limits = c(0, 50)
    ) +
    scale_x_discrete(
      breaks = c("cv", "cv_inter", "cv_intra"),
      labels = c("total", "between", "within")
    ) +
    theme(
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = 0.5,
        face = "plain"
      ),
      axis.text.y = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = .5,
        face = "plain"
      ),
      #plot.margin = margin(2, 2, 2, 2, "cm"),
      strip.text.x = element_text(size = 12),
      #legend.position=c(.10, .98),
      #legend.title = element_blank(),
      legend.position = "right",
      axis.title.y = element_text(
        size = rel(1.3),
        angle = 0,
        vjust = 0.5
      ),
      axis.title.x = element_text(size = rel(1.3), angle = 0),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      # switch off major gridlines
      panel.grid.minor = element_blank(),
      # switch off minor gridlines
      panel.border = element_blank(),
      axis.line = element_line(size = 0.5),
      #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
      #panel.border = element_border(c("left","bottom")),
      legend.title = element_blank(),
      # switch off the legend title
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.0, "lines")
      
    )
  
  for (x in colnames(t_esp.i)) {
    #x=colnames(t_esp.i)[13]
    is_na = (is.na(t_esp.i[, x]) | t_esp.i[, x] == 0)
    if (sum(is_na) > 0) {
      na_rownames = gsub("\\..*", "", rownames(t_esp.i[is_na,]))
      na_sample_group = gsub("(\\D+)(\\d+$)", "\\1", na_rownames)
      for (y in unique(na_sample_group)) {
        #y=na_sample[1]
        t_esp.i[which(is_na)[grepl(y, na_sample_group)], x] =
          rnorm(
            n = sum(grepl(y, na_sample_group)),
            mean = median(t_esp.i[grepl(y, rownames(t_esp.i)), x], na.rm =
                            T),
            sd = median(t_esp.i[grepl(y, rownames(t_esp.i)), x], na.rm =
                          T) * 0.01
          )
      }
    }
    
  }
  # Analisis PCA especies
  
  t_QC.esp.sc <- scale(t_esp.i, scale = T)
  
  t_QC.esp.svd <- svd(t_QC.esp.sc)
  t_QC.esp.scores <- t_QC.esp.svd$u %*% diag(t_QC.esp.svd$d)
  t_QC.esp.loadings <- t_QC.esp.svd$v
  t_QC.esp.vars <- t_QC.esp.svd$d ^ 2 / (nrow(t_QC.esp.sc) - 1)
  t_QC.esp.totalvar <- sum(t_QC.esp.vars)
  t_QC.esp.relvars <- t_QC.esp.vars / t_QC.esp.totalvar
  variances <- 100 * round(t_QC.esp.relvars, digits = 3)
  variances[1:5]
  #sum(variances[1:2])
  components <- c(1, 2)
  ####GRAFICOS PCA DE las tandas #####
  ###################################
  par(
    mfrow = c(1, 2),
    mar = c(4, 4, 2, 1),
    oma = c(1, 0.0, 1, 0),
    cex = 1.2
  )
  
  
  #Verificar que no se han eliminado condiciones
  t_esp.i <-
    cbind(t_esp[rownames(t_esp.i), c("type", "tanda")], t_esp.i)
  dfxy <-
    cbind(
      type = t_esp.i$type,
      tanda = t_esp.i$tanda,
      data.frame(t_QC.esp.scores[, 1:2])
    )
  # save(dfxy,escala_batch,file="qc_plot.rda")
  # load("qc_plot.rda")
  #
  
  dfxy$label = NA
  labelToShow = match(unique(dfxy$type), dfxy$type)
  dfxy[labelToShow, "label"] = as.character(unique(dfxy$type))
  dfxy$opacity = ifelse(dfxy$type == "QC", 1, 0.9)
  dfxy$xaxis = ifelse(dfxy$type == "QC", min(dfxy$X2) * 0.9, median(dfxy$X2) *
                        1.3)
  #View(dfxy)
  graf_p1 <- ggplot(dfxy, aes(x = X2, y = X1)) +
    geom_hline(yintercept = 0,
               color = "grey",
               linetype = 'dashed') +
    geom_vline(xintercept = 0,
               color = "grey",
               linetype = 'dashed') +
    geom_point(
      aes(
        shape = type,
        fill = as.character(tanda),
        colour = as.character(tanda),
        alpha = opacity
      ),
      #colour="black",
      size = 3,
      stroke = 0.75,
      #alpha = .8,
      show.legend = T
    ) +
    
    scale_shape_manual(values = c(21:25)) +
    
    scale_fill_manual(values = escala_batch) +
    scale_color_manual(values = escala_batch) +
    stat_ellipse(
      aes(x = X2, y = X1, shape = type),
      color = "darkgrey",
      alpha = 1,
      show.legend = F
    ) +
    labs(
      title = "PCA QC vs SAMPLES",
      x = paste("PC ", components[2], " (", variances[components[2]], "%)", sep = ""),
      y = paste("PC ", components[1], " (", variances[components[1]], "%)", sep = "")
    ) +
    geom_label(
      #label.size = 0,
      position = position_dodge(0.5),
      aes(label = label,
          x = xaxis),
      #label.padding = unit(0.1, "lines"),
      #label.r	= unit(1, "lines"),
      color = "black",
      size = 5,
      show.legend = FALSE
    ) +
    geom_point(
      aes(
        shape = type,
        fill = as.character(tanda),
        colour = as.character(tanda),
        alpha = opacity
      ),
      colour = "black",
      size = 3,
      stroke = 0.75,
      #  alpha = .8,
      show.legend = F
    ) +
    # labs(
    #   title = "PCA QC vs SAMPLES",
    #   x = paste("PC ", components[2], " (", variances[components[2]], "%)", sep = ""),
    #   y = paste("PC ", components[1], " (", variances[components[1]], "%)", sep = "")
    # ) +
    guides(
      colour = guide_legend(override.aes = list(colour = escala_batch)),
      shape =
        "none",
      alpha = "none"
    ) +
    theme(
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = 0.5,
        face = "plain"
      ),
      axis.text.y = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = .5,
        face = "plain"
      ),
      strip.text.x = element_text(size = 12),
      legend.position = "top",
      axis.title.y = element_text(
        size = rel(1.3),
        angle = 90,
        vjust = 0.0
      ),
      axis.title.x = element_text(size = rel(1.3), angle = 0),
      strip.background = element_blank(),
      axis.line = element_line(size = 0.5),
      legend.title = element_blank(),
      # switch off the legend title
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.0, "lines"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      )
    )
  
  # graf_p1
  # resid_1
  dfxy$tanda <- as.character(dfxy$tanda)
  dfxy[which(dfxy$type == "QC"), "tanda"] <- c("QC")
  #dfxy$subgroups4<-factor(dfxy$subgroups4,levels=as.character(sort(unique(dfxy$subgroups4),decreasing = T)),ordered=T)
  dfxy$tanda <-
    factor(dfxy$tanda,
           levels = c("QC", as.character(c(1:20))),
           ordered = T)
  
  graf_p2 <- ggplot(dfxy, aes(x = X2, y = X1)) +
    geom_hline(yintercept = 0,
               color = "grey",
               linetype = 'dashed') +
    geom_vline(xintercept = 0,
               color = "grey",
               linetype = 'dashed') +
    geom_point(
      aes(shape = tanda, fill = tanda),
      size = 3,
      stroke = 0.75,
      alpha = .8
    ) +
    scale_shape_manual(values = c(21:25)) +
    scale_fill_manual(values = escala) +
    stat_ellipse(aes(x = X2, y = X1, shape = tanda),
                 color = "darkgrey",
                 alpha = .8) +
    #stat_ellipse(aes(x=X1,y=X2,shape=subgroups[2]))+
    labs(
      title = "Between batch dispersion",
      x = paste("PC ", components[2], " (", variances[components[2]], "%)", sep = ""),
      y = paste("PC ", components[1], " (", variances[components[1]], "%)", sep = "")
    ) +
    theme(
      axis.text.x = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = 0.5,
        face = "plain"
      ),
      axis.text.y = element_text(
        size = 12,
        angle = 0,
        hjust = .5,
        vjust = .5,
        face = "plain"
      ),
      strip.text.x = element_text(size = 12),
      legend.position = "top",
      axis.title.y = element_text(
        size = rel(1.3),
        angle = 90,
        vjust = 0.0
      ),
      axis.title.x = element_text(size = rel(1.3), angle = 0),
      strip.background = element_blank(),
      axis.line = element_line(size = 0.5),
      legend.title = element_blank(),
      # switch off the legend title
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.0, "lines"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      )
    )
  
  QC_data <- list(
    dfc.qc.es.total,
    #1
    dfc.qc.sum.total,
    #2
    t_SUM.i,
    #3
    graf_cv_especies,
    #4
    graf_cv_clases,
    #5
    graf_p1,
    #6
    graf_p2,
    #7
    resid_1 #8
  )
  names(QC_data) <- c(
    "cv_especies",
    "cv_sumas",
    "t_SUM_imputadas",
    "graf_cv_especies",
    "graf_cv_clases",
    "PCA_QC",
    "PCA_QC2",
    "graf_residuos"
  )
  return(QC_data)
}


PCA_analysis<-function(t_esp.i,casos,medidas,groupVar,rotate=F,invert=F){
  ##recibe una tablas y calcula la calidad de la cuantificacion
  ## La tabla debe incluir SAMPLEid = (identificacion de la muestra)
  ##                       tanda = si hay varias tandas de analisis
  ##                       type  = si es control (QC) o muestra (SAMPLE)
  ## Las imputaciones se realizan seg?n rownames que tiene que tener la estructura (GRUPO_sAMPLEid)
  
  #t_esp<-t_rsALL.concTOistd
  require('mixOmics')
  require(ade4)
  #t_esp.i<-t_rsRESUM.sum
  #eindex<-rsALL.concTOistd[,c("index","NAME_ltr","avg_mz","avg_Sec","Region")]
  #casos<-cases
  #medidas<-measures.sum
  #groupVar<-"Fsig"
  #rotate=F
  #invert=F
  
  t_esp.i<-t_esp.i[which(t_esp.i$type %in% c("SAMPLE")),c(casos,medidas)]
  #t_esp.QC<-cbind(SAMPLEid=rownames(t_esp.QC),t_esp.QC)
  
  # Analisis PCA 
  t_QC.esp.sc<-t_esp.i[,which(colnames(t_esp.i) %in% medidas)]#seleccionamos medidas
  t_QC.esp.sc<-scale(t_QC.esp.sc,scale=T);
  t_QC.esp.svd <- svd(t_QC.esp.sc)
  t_QC.esp.scores <- t_QC.esp.svd$u %*% diag(t_QC.esp.svd$d)
  t_QC.esp.loadings <- t_QC.esp.svd$v
  t_QC.esp.vars <- t_QC.esp.svd$d^2/(nrow(t_QC.esp.sc) - 1)
  t_QC.esp.totalvar <- sum(t_QC.esp.vars)
  t_QC.esp.relvars <- t_QC.esp.vars / t_QC.esp.totalvar
  variances <- 100 * round(t_QC.esp.relvars, digits = 3)
  variances[1:5]
  sum(variances[1:2])
  components<-c(1,2)
  ####GRAFICOS PCA DE las tandas #####
  ###################################
  escala<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#fee090','#ffff33','#a65628','#f781bf')
  escala<-escala[c(1,5,2,7,3,6,4,8)]  
  dfxy<-cbind(group=t_esp.i[,which(colnames(t_esp.i) %in% groupVar)],data.frame(t_QC.esp.scores[,1:5]))
  rownames(dfxy)<-rownames(t_esp.i)
  if((rotate)==T){colnames(dfxy)[c(2,3)]<-colnames(dfxy)[c(3,2)]}
  if((invert)==T){dfxy$X1<-dfxy$X1*c(-1)}
  
  graf_PCA<-ggplot(dfxy, aes(x=X1, y=X2)) +
    geom_hline(yintercept=0,color="grey",linetype = 'dashed')+
    geom_vline(xintercept=0,color="grey",linetype = 'dashed')+
    geom_point(aes(shape=group,fill=group),size=3,stroke=0.75,alpha=.8)+
    scale_shape_manual(values=c(21:25))+
    scale_fill_manual(values=escala)+
    stat_ellipse(aes(x=X1,y=X2,shape=group),color="darkgrey",alpha=.8)+
    labs(title="PCA plot", 
         x = paste("PC ",components[1]," (",variances[components[1]], "%)", sep = ""),
         y = paste("PC ",components[2]," (",variances[components[2]], "%)", sep = ""))+
    theme(axis.text.x=element_text(size=12,angle=0,hjust=.5,vjust=0.5,face="plain"),
          axis.text.y=element_text(size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          strip.text.x=element_text(size=12),
          legend.position="top",
          axis.title.y = element_text(size = rel(1.3), angle = 90, vjust = 0.0),
          axis.title.x = element_text(size = rel(1.3), angle = 0),
          strip.background=element_blank(),
          axis.line = element_line(size=0.5),
          legend.title = element_blank(), # switch off the legend title
          legend.text = element_text(size=12),
          legend.key.size = unit(1.0, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )
  
  graf_PCA
  
  ###GRAFICO DE LOADINGS
  
  
  rownames(t_QC.esp.loadings)<-medidas
  
  #rownames(t_QC.esp.loadings)<-eindex[which(eindex$index %in% medidas),"NAME_ltr"]
  
  #eindex[which(eindex$index %in% medidas),]
  #loads<-cbind(eindex[which(eindex$index %in% medidas),],data.frame(t_QC.esp.loadings))
  #loads<-cbind(group=t_esp.i[,which(colnames(t_esp.i) %in% groupVar)],data.frame(t_QC.esp.scores[,1:2]))
  
  t_QC.esp.loadings<-data.frame(t_QC.esp.loadings)
  t_QC.esp.loadings<-cbind(Region=gsub("(.*)\\_(.*)","\\1",rownames(t_QC.esp.loadings)),t_QC.esp.loadings)
  
  #loadings$tanda<-as.character(dfxy$group)
  #dfxy[which(dfxy$type=="QC"),"tanda"]<-c("QC")
  #dfxy$subgroups4<-factor(dfxy$subgroups4,levels=as.character(sort(unique(dfxy$subgroups4),decreasing = T)),ordered=T)
  #dfxy$tanda<-factor(dfxy$tanda,levels=c("QC",as.character(c(1:20))),ordered=T)
  #escala<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#fee090','#ffff33','#a65628','#f781bf')
  #'#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'
  escala2<-c(brewer.pal(12,"Paired"),"black")
  
  graf_LOADS<-ggplot(t_QC.esp.loadings, aes(x=X1, y=X2,fill=Region)) +
    geom_hline(yintercept=0,color="grey",linetype = 'dashed')+
    geom_vline(xintercept=0,color="grey",linetype = 'dashed')+
    geom_segment(aes(x = 0, y = 0, xend = X1, yend = X2,colour=Region),
                 lineend = "round", # See available arrow types in example above
                 linejoin = "round",
                 size = 0.5, 
                 arrow = arrow(length = unit(0.2, "cm"))
    )+
    geom_text(aes(color=Region),nudge_x=0.01,nudge_y=0.01, 
              check_overlap=T,
              label=rownames(t_QC.esp.loadings),size=3, show.legend = F)+
    #geom_point(shape=21,size=3,stroke=0.75,alpha=.8)+
    scale_fill_manual(values=escala2)+
    scale_color_manual(values=escala2)+
    #stat_ellipse(aes(x=X2,y=X1,shape=tanda),color="darkgrey",alpha=.8)+
    #stat_ellipse(aes(x=X1,y=X2,shape=subgroups[2]))+
    labs(title="Loadings plot", 
         x = paste("PC ",components[1]," (",variances[components[1]], "%)", sep = ""),
         y = paste("PC ",components[2]," (",variances[components[2]], "%)", sep = ""))+
    theme(axis.text.x=element_text(size=12,angle=0,hjust=.5,vjust=0.5,face="plain"),
          axis.text.y=element_text(size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          strip.text.x=element_text(size=12),
          legend.position="top",
          axis.title.y = element_text(size = rel(1.3), angle = 90, vjust = 0.0),
          axis.title.x = element_text(size = rel(1.3), angle = 0),
          strip.background=element_blank(),
          axis.line = element_line(size=0.5),
          legend.title = element_blank(), # switch off the legend title
          legend.text = element_text(size=12),
          legend.key.size = unit(1.0, "lines"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )
  
  graf_LOADS
  
  #### EXPLORAR PLS-DA ###
  
  
  X<-t_QC.esp.sc
  dim(X)
  Y<-t_esp.i[,which(colnames(t_esp.i) %in% groupVar)]
  
  length(Y)==dim(X)[1]
  
  MyResult.plsda<-plsda(X,Y,ncomp=3)
  graf_PLSDA<-plotIndiv(MyResult.plsda, ind.names = F, legend=TRUE, comp=c(1,2), pch=c(16),cex=3,style = 'ggplot2',
                        ellipse = TRUE, star = F, title = 'PLS-DA',col.per.group=escala[1:length(levels(Y))],
                        X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')
  graf_PLSDA
  graf_PLSDA_vars<-plotVar(MyResult.plsda, var.names=T, col ="black", comp = c(1,2))
  graf_PLSDA_vars
  auc.plsda <- auroc(MyResult.plsda,roc.comp=1)
  auc.plsda
  #auc.plsda <- auroc(MyResult.plsda2)
  
  PCA_data<-list(dfxy,t_QC.esp.loadings,graf_PCA,graf_LOADS,graf_PLSDA,auc.plsda)
  names(PCA_data)<-c("pca_scores","pca_loadings","graf_PCA","graf_loadings","graf_plsda","graf_auc")
  return(PCA_data)
  
  
}

#### graficas pvalue
summaryRES_PACK_2<-function(dfc_obj, digitos=2, varTOcast, lipidclassColName=NULL, estat="media_se",my_comparisons,salto=1.1){
  
  #funcion recibe un objeto dfc y requiere columna "variable" para las classes
  #digitos=2 numero de digitos de redondeo
  #varTOcast nombre del la variable que agrupa los casos segun tratamiento
  #salto es el multiplicador de la altura maxima de la grafica
  #estat, por defecto la media_sd, mediana_ci estat="mediana", media_se, estat=media_se
  #lipidclassColName para cuando es distinto a "variable"
  
  ## test values
  #dfc_obj<-dfc.fibro;
  #digitos<-2;
  #dfc_obj<-dfc_AUC
  #my_comparisons<-mapply(cbind,c("CNT","HFD","MCD"),c("CNT_FEN","HFD_FEN","MCD_FEN"),SIMPLIFY = F)
  
  # dfc_obj=dfc.sum.sel
  # digitos=2
  # varTOcast = "GRUPO"
  # estat="media_se"
  # my_comparisons=comparaciones
  # browser()
  dfc_obj$msd<-paste0(paste(round(dfc_obj$value,digits=digitos),round(dfc_obj$sd,digits=digitos),sep=" \u00b1 ")," (",dfc_obj$N,")")
  dfc_obj$medci<-paste0(paste(round(dfc_obj$median,digits=digitos),round(dfc_obj$ci,digits=digitos), sep=" \u00b1 ")," (",dfc_obj$N,")")
  dfc_obj$mse<-paste0(paste(round(dfc_obj$value,digits=digitos),round(dfc_obj$se,digits=digitos), sep=" \u00b1 ")," (",dfc_obj$N,")")
  
  # guardar todas las variables del entorno, para arreglar errores
  # save(list = ls(),
  #      file = "summaryRES_PACK_2.Rdata",
  #      envir = environment())
  # load("summaryRES_PACK_2.Rdata")
  browser()
  if (is.na(my_comparisons)) {
    my_groups<-unique(dfc_obj[,1])
    dfx<-expand.grid(my_groups,my_groups,stringsAsFactors = T)
    dfx<-dfx[as.character(dfx$Var1) <= as.character(dfx$Var2),]
    dfx<-dfx[!dfx$Var1==dfx$Var2,]
    dfx<-dfx[order(dfx[,1]),]
    my_comparisons<-mapply(cbind, as.character(dfx$Var1), as.character(dfx$Var2), SIMPLIFY = F)
    
  }else {
    my_comparisons<-my_comparisons
  }
  
  p<-NULL
  p<-list();  
  p.y.position=NULL
  p.y.position=list()
  ## Comparaciones por pares 
  for (i in 1:length(my_comparisons)){  
    browser()
    #i=1
    MEAN1<-dfc_obj[which(dfc_obj[,varTOcast] %in% my_comparisons[[i]][,1]),"value"]  
    MEAN2<-dfc_obj[which(dfc_obj[,varTOcast] %in% my_comparisons[[i]][,2]),"value"]
    SD1<-dfc_obj[which(dfc_obj[,varTOcast] %in% my_comparisons[[i]][,1]),"sd"]  
    SD2<-dfc_obj[which(dfc_obj[,varTOcast] %in% my_comparisons[[i]][,2]),"sd"]
    N1<-dfc_obj[which(dfc_obj[,varTOcast] %in% my_comparisons[[i]][,1]),"N"]  
    N2<-dfc_obj[which(dfc_obj[,varTOcast] %in% my_comparisons[[i]][,2]),"N"]  
    tst<-(MEAN1-MEAN2)/sqrt((SD1^2/N1)+ (SD2^2/N2))
    df<-((SD1^2/N1)+(SD2^2/N2))^2/(((SD1^2/N1)^2/(N1-1))+((SD2^2/N2)^2/(N2-1)))
    p[[i]]<-2*pt(-abs(tst), df) #corregido 
    names(p)[i]<- paste0(my_comparisons[[i]][,1],"_vs_",my_comparisons[[i]][,2])
    
    if(!is.null(lipidclassColName)){
      colnames(dfc_obj)[which(colnames(dfc_obj) %in% lipidclassColName)]="variable"
    }
    
    ###creamos una lista con valores maximos
    sub_dfc_obj=dfc_obj
    sub_dfc_obj[is.nan(sub_dfc_obj$value),c('value','median','cv','se','ci')]=0
    sub_dfc_obj[which(sub_dfc_obj[,varTOcast] %in% my_comparisons[[i]]),varTOcast]=sub_dfc_obj[which(sub_dfc_obj[,varTOcast] %in% my_comparisons[[i]]),varTOcast]
    get_max_yscale=aggregate(as.formula(paste0(".~variable+",varTOcast)),sub_dfc_obj[which(sub_dfc_obj[,varTOcast] %in% my_comparisons[[i]]),c("variable","value","sd",varTOcast)],sum)
    
    #get_max_yscale=aggregate(.~variable,sub_dfc_obj[,c("variable","value","sd")],sum)
    p.y.position[[i]]=get_max_yscale$value+get_max_yscale$sd  ### para evitar que se sobreponga a las barras de error
    names(p.y.position)[i]<- paste0(my_comparisons[[i]][,1],"_vs_",my_comparisons[[i]][,2])
  }  
  
  if (estat=="mediana") {
    rsdfc_obj<-reshape2::dcast(dfc_obj, paste("variable", "~",varTOcast),value.var="medci") 
  } else if ( estat=="media_se") {
    rsdfc_obj<-reshape2::dcast(dfc_obj, paste("variable", "~",varTOcast),value.var="mse")
  } else {
    rsdfc_obj<-reshape2::dcast(dfc_obj, paste("variable", "~",varTOcast),value.var="msd")
  }
  
  kk<- t(data.frame(t(sapply(p, unlist))))
  p_y.position<<-t(data.frame(t(sapply(p.y.position, unlist))))
  
  rsdfc_obj<-cbind(rsdfc_obj,kk)
  
  #check<-rsdfc_obj[,c(1:2,4,8)]
  
  #### creamos una variable global para poner p.signf
  get_max_yscale=aggregate(.~variable,dfc_obj[,c("variable","value","sd")],max)
  
  p_y.position=cbind("V1"=as.character(rsdfc_obj[,1]),as.data.frame(p_y.position,stringsAsFactors = F))
  p_y.position[,2]=get_max_yscale$value+(get_max_yscale$sd*2)
  if (ncol(p_y.position)>2) { #a?adido para evitar error en comparaciones de solo dos grupos 
  for(i in c(3:ncol(p_y.position))){
    #vamos sumando un porcentaje de altura maxima a cada nueva bracket de pvalues
    #i=2
    p_y.position[,i]=p_y.position[,(i-1)]+((salto-1)*p_y.position[,2])
  }
  }else{
    p_y.position[,2]=p_y.position[,2]+p_y.position[,2]*(salto-1)
  }
  p_y.position=reshape2::melt(p_y.position)
  
  kk<-as.data.frame(cbind(as.character(rsdfc_obj[,1]),kk),stringsAsFactors=FALSE)
  p_table<-reshape2::melt(kk,id.vars="V1")
  if(length(unique(p_table$V1))<2){
    p_table$variable=rownames(kk)
  }
  p_table$group1<-sub("(.*)(_vs_)(.*)","\\1",p_table[,2])
  p_table$group2<-sub("(.*)(_vs_)(.*)","\\3",p_table[,2])
  #colnames(p_y.position)=colnames(p_table)[c(1,3)]
  p_table<-unique(merge(p_table,p_y.position,by=c("V1","variable")))
  p_table$value.x<-as.numeric(p_table$value.x) #Fix error
  p_table$p.signif<-ifelse(is.na(p_table$value.x),"ns",
                           ifelse(p_table$value.x<=0.0001,"****",
                                  ifelse(p_table$value.x<=0.001,"***",
                                         ifelse(p_table$value.x<=0.01,"**",
                                                ifelse(p_table$value.x<0.05,"*",
                                                       ifelse(p_table$value.x>0.1,"ns",p_table$value.x)
                                                )))))
  colnames(p_table)<-c("supp",".y.","p.format","group1","group2","y.position","p.signif")
  p_table$.y.<-"value"
  #p_table<<-p_table
  #print("Se ha creado una variable p_table con datos necesarios para indicar p-value en las graficas")
  ####
  #list(stat=rsdfc_obj,p_table=p_table)
  return(list(stat=rsdfc_obj,p_table=p_table)) 
  
  
}

