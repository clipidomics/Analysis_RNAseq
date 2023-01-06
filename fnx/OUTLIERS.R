####outliers bohdan
cv_outlier <- function(info_dataframe,r=30) {
  ##devuelve un vector de T/F
  cv = sd(info_dataframe, na.rm = TRUE)/median(info_dataframe, na.rm = TRUE)*100
  mediana_inicio=median(info_dataframe, na.rm = TRUE)
  not_NAss=length(info_dataframe)-sum(is.na(info_dataframe))
  if (!is.na(cv) & cv > r & not_NAss>2) {
    #if(exists("outlier")) rm(outlier)
    while (cv > r & not_NAss>2) {
      
      if(exists("outlier")){
        outlier_resta = abs(info_dataframe[which(outlier!=T)]-mediana_inicio)
        #corregimos NAs despues de la resta para que no sean seleccionados como outlier
        outlier_resta[is.na(outlier_resta)]=0
        outlier[which(outlier!=T)] = as.vector(outlier_resta==max(outlier_resta))
      } else {
        outlier_resta = abs(info_dataframe-mediana_inicio)
        #corregimos NAs despues de la resta para que no sean seleccionados como outlier
        outlier_resta[is.na(outlier_resta)]=0
        outlier = as.vector(outlier_resta==max(outlier_resta))
      }
      
      
      cv = sd(info_dataframe[which(outlier!=T)], na.rm = TRUE)/mediana_inicio*100
      not_NAss=length(info_dataframe[which(outlier!=T)])-sum(is.na(info_dataframe[which(outlier!=T)]))
    }
    
  }   else  {
    info_dataframe[1:length(info_dataframe)]=F
    outlier=info_dataframe
  }
  return(outlier)
}

new_outlier <- function(info_dataframe,
                        Rval=2,
                        digits=2,
                        method="hb",
                        ISTD_variabily=ISTD_variabilidad,
                        index=rsALL[,c("mz_rt","NAME")],
                        imput_max=F,
                        umbral=10,
                        informe=F,
                        imput_mean=T,
                        sel_columns=F,
                        sel_rows=F,
                        regex_filas="_[0-9]*$",
                        umbral_factor=40,
                        fig_col=4,
                        fig_row=4,
                        dotsize=3,
                        desviacion_ISTD=15,
                        porcentaje_resumen=5,
                        NAs_for_median=F,
                        custom_imput_sd=F,
                        use_mean=F){
  # # # # FORMATO DE LOS PARAMETROS
  # info_dataframe=t_rsALL.concTOcal[,-c(1:2)] #matriz de filas para grupos y columnas para especies
  # Rval=3 #2 para hb lo que significa 2 sd o un numero sin porcentaje para CV
  # digits=2 #digitos para los numeros que salen en los mensajes
  # method="hb" #metodo "hb" para doble sd o "CV" para el coeficiente de variacion
  # umbral=10 #umbral para determinar en que porcentaje de especies tiene que imputarse para considerarse una muestra extrema.
  # informe=F #imprimir un informe con cada muestra que ha sido imputada. T o F
  # imput_mean=T #si imputar "T" o eliminar "F" los outliers. influye en las sumas
  # sel_columns=F #si queremos imputas solo algunas columnas ponemos el numero "1:10" para seleccionar una submatriz o FALSE para desactivar
  # sel_rows=F #si queremos imputas solo algunas filas ponemos el numero "1:10" para seleccionar una submatriz o FALSE para desactivar
  # regex_filas="_.*$" #regex para hacer grupos
  # fig_col=4 #columnas para las graficas de outliers imputados
  # fig_row=4 #filas para lo mismo
  # dotsize=3 #tamaÃ±o de los puntos para lo mismo
  # imput_max=F ###numero de valores que quieres imputar porcentaje de muestras que se permiten imputar, enter 0 y 1
  # ISTD_variabily=F ### un dataframe de valores de los ISTD con nombre de filas iguales al los nombres de ttos en info_dataframe o false cuando no lo haya
  # index=eindex[,c("index","NAME_ltr")] ###coge dos columnas, la primera con nombre de columna de info_dataframe y la segunda su correspondencia para pintar las graficas
  # porcentaje_resumen=5 #para el informe
  # desviacion_ISTD=20 #para el informe
  # NAs_for_median=F #si quieres sustituir todos los NAs por la mediana de su grupo
  # custom_imput_sd=0.01 # es el valor por el cual quieres multiplicar la mediana para disminuir variabilidad durante la imputacion
  # use_mean=F
  ### rm(n_imput_max,same_group,index,info_dataframe,Rval,digits,method,umbral,informe,imput_mean,sel_columns,sel_rows,regex_filas)
  ##los datos tienen que estar agrupados por filas
  log_file=paste0("log_imputaciones_",as.numeric(Sys.time()),".txt")
  sink(file=log_file,split = T)
  print(Sys.time())
  print(getwd())
  
  print(paste("Rval",Rval,
              "method",method,
              "imput_max",imput_max,
              "umbral",umbral,
              "imput_mean",imput_mean,
              "sel_columns",sel_columns,
              "sel_rows",sel_rows,
              "regex_filas",regex_filas,
              "umbral_factor",umbral_factor,
              "fig_col",fig_col,
              "fig_row",fig_row,
              "dotsize",dotsize,
              "porcentaje_resumen",porcentaje_resumen,
              "NAs_for_median",NAs_for_median,
              "custom_imput_sd",custom_imput_sd,
              "use_mean",use_mean,sep = " | "))
  #hacemos una copia por si hay que quitar columnas o filas
  original_info_dataframe=info_dataframe
  #filtramos columnas/filas
  if(sum(sel_columns!=F)>0 & sum(sel_rows!=F)>0){
    info_dataframe=info_dataframe[sel_rows,sel_columns]
  } else if(sum(sel_columns!=F)>0 & sum(sel_rows!=F)==0){
    info_dataframe=info_dataframe[,sel_columns]
  } else if(sum(sel_columns!=F)==0 & sum(sel_rows!=F)>0){
    info_dataframe=info_dataframe[sel_rows,]
  }
  #hacemos una copia de la matriz de los datos que se utilizara en caso de que no queramos sustituir todos los outliers
  info_dataframe_backup<<-info_dataframe
  
  ##agrupamos por filas
  row_names=gsub(regex_filas,"",rownames(info_dataframe))
  
  
  ######    IMPUTACION    #########
  
  ##vamos de fila en fila
  for(rows in unique(row_names)) {
    #rows=unique(row_names)[2]
    #rows="ND"
    ##como minimo tiene que haber 3 repeticiones para realizar la busquedade los outliers
    if((table(row_names)[rows])>2) {
      
      ##buscamos outliers columna por columna
      for(col in colnames(info_dataframe)){
        #col=colnames(info_dataframe)[1]
        #
        ###si queremos quitar todos los NAs
        if(NAs_for_median==T){
          info_dataframe[row_names %in% rows,col][(info_dataframe[row_names %in% rows,col]==0)]=NA
          if(custom_imput_sd!=F){
            imput_mean_sd=median(info_dataframe[row_names %in% rows,col],na.rm=T)*custom_imput_sd
          }else{
            imput_mean_sd=median(info_dataframe[row_names %in% rows,col],na.rm=T)*0.05
          }
          for (k in c(which(is.na(info_dataframe[row_names %in% rows,col]) ))){
            info_dataframe[row_names %in% rows,col][k]=rnorm(1, mean=median(info_dataframe[row_names %in% rows,col],na.rm=T), sd=imput_mean_sd)
          }
        }
        
        ##hacemos una copia para el informe
        ##si hay mas de 3 datos no NAs
        if(length(info_dataframe[row_names %in% rows,col][!is.na(info_dataframe[row_names %in% rows,col])])>2){
          
          previous=format(info_dataframe_backup[row_names %in% rows,
                                                col],digits = digits)
          
          
          ##buscamos los outliers, seleccionamos solo el mas alejado de la mediana
          if(tolower(method)=="cv") {
            #example=matrix(info_dataframe[row_names %in% rows,col])
            outlier=cv_outlier(matrix(info_dataframe[row_names %in% rows,col]),r=Rval)
          } else {
            outlier=hboutlier(matrix(info_dataframe[row_names %in% rows,col]),r=Rval)
          }
          
          
          
          
          ## si hay outlier
          if(sum(outlier)>0){
            #seleccionamos el mas alejado
            ###sustituimos los NAs por un valor muy muy alejado de la mediana para que sean primeros en sustituirse
            if(imput_max!=F){
              n_imputs_allowed=round(sum(grepl(rows,row_names))*imput_max)
            } else {
              n_imputs_allowed=round(sum(grepl(rows,row_names))*0.30)
            }
            NAs_detectados=rep(NA,length(is.na(info_dataframe[row_names %in% rows,col][outlier]))) ## creamos un vector de NAs detectados para luego sustituir los que se permiten imputar segun el numero de casos por un valor muy alto para que sean detectados entre primeros como oulier
            if(length(NAs_detectados)<n_imputs_allowed){
              n_imputs_allowed=length(NAs_detectados)
            }
            
            ###solo trabajamos con los valores, no tocamos los NAs
            outlier[is.na(info_dataframe[row_names %in% rows,col])]=FALSE
            
            ##buscamos los mas alejados de la mediana
            selected_outlier=info_dataframe[row_names %in% rows,col][outlier]-median(info_dataframe[row_names %in% rows,col],na.rm=T)
            selected_outlier=match(sort(abs(selected_outlier),decreasing = T),abs(selected_outlier))
            
            ##filtramos para seleccionar el mas alejado de la mediana si hay pocas repeticiones
            
            outlier_filtro=outlier
            outlier_filtro[1:length(outlier_filtro)]=F
            outlier_filtro[which(outlier==T)[selected_outlier[1:n_imputs_allowed]]]=T
            outlier=outlier_filtro
            
            
            
            #calculamos la media para imputar excluyendo los outliers y las NAs
            if(use_mean==T){
              imput_mean_value=mean(matrix(info_dataframe[row_names %in% rows,col])[-which(outlier)],na.rm = T)
              imput_mean_sd=sd(matrix(info_dataframe[row_names %in% rows,col])[-which(outlier)],na.rm = T)
            } else {
              ### comprobamos, si los datos son "normales" usamos la media, sino la mediana
             
                imput_mean_value=median(matrix(info_dataframe[row_names %in% rows,col]),na.rm = T)
                imput_mean_sd=sd(matrix(info_dataframe[row_names %in% rows,col])[-which(outlier)],na.rm = T)
                
              
              
            }
            
            ### imputaciones
            if(custom_imput_sd!=F){
              imput_mean_sd=imput_mean_value*custom_imput_sd
            }
            
            if(informe==T){
              ##imprimimos el informe
              print(cbind(Columna=col,
                          Fila=rownames(info_dataframe)[row_names %in% rows],
                          Valor=previous,
                          Outlier=outlier,deparse.level = 0,Media=format(imput_mean_value,digits = digits),SD=format(imput_mean_sd,digits = digits)))
              print(paste0("--------------------------------------"))
            }
            
            for(h in which(outlier==T)){
              if(imput_mean==T){
                info_dataframe[row_names %in% rows,col][h]=rnorm(1, mean=imput_mean_value, sd=imput_mean_sd)
              } else {
                info_dataframe[row_names %in% rows,col][h]=NA
              }
            }
          } 
        }
      }
    }
  }
  
  
  ######    GRAFICAS  e INFORME  #########
  
  
  
  #### PREPARAMOS DATOS PARA HACER LAS GRAFICAS
  if(is.data.frame(index)==T){
    #preparamos los datos en long format hacer la grafica
    plot_data=reshape2::melt(cbind("ttos"=rownames(info_dataframe_backup),info_dataframe_backup))
    plot_data_imputado=reshape2::melt(cbind("ttos"=rownames(info_dataframe),info_dataframe))
    plot_data=cbind(plot_data,"grupo"=gsub(paste0("(.*)",regex_filas),"\\1",plot_data$ttos),"caso"=gsub("(.*)_(.*$)","\\2",plot_data$ttos),"value_imput"=plot_data_imputado$value,"value_old"=plot_data$value)
    
    ##creamos una variable con todas las imputaciones
    plot_data_imputado=which(plot_data$value_imput!=plot_data$value)
    NAs_antes=which(is.na(plot_data$value_imput)!=is.na(plot_data$value))###cogemos las NAs imputadas que no coinciden con los que habia nates
    plot_data_imputado=append(plot_data_imputado,NAs_antes)
    plot_data_imputado=unique(plot_data_imputado)
    
    plot_data[!(c(1:nrow(plot_data)) %in% plot_data_imputado),"value_imput"]=NA ### quitamos valores no sustituidos
    plot_data[!(c(1:nrow(plot_data)) %in% plot_data_imputado),"value_old"]=NA ### quitamos valores antiguos no sustituidos
    plot_data[(c(1:nrow(plot_data)) %in% plot_data_imputado),"value"]=NA ### quitar valores sustituidos
    
    
    
    #casamos mzrt con la region correspondiente
    regions_id=(cbind("variable"=as.character(index[,1]),"NAME"=as.character(gsub("(.*)#.*","\\1",index[,2]))))
    #regions_id[] <- lapply(regions_id, as.character)
    plot_data=merge(plot_data,regions_id)
    classes=unique(plot_data[plot_data_imputado,"NAME"])
    plot_data=plot_data[(plot_data$NAME %in% classes),]
    plot_data=plot_data[order(plot_data$NAME),]
    
  }
  ###grafica movida hasta despues de impirmir los resultados
  if(!exists("remove_rows")) remove_rows=0
  
  print(paste0("Resumen de condiciones con mas de ",porcentaje_resumen,"% de especies outlier:"))
  outliers_resumen_sustituciones  <<-data.frame()
  
  ###para casos cuando no tenemos tabla de istds
  if(ISTD_variabily==F){
    info_dataframe_backup_fictiona_istds=info_dataframe_backup
    colnames(info_dataframe_backup_fictiona_istds)=sub("(.*)_.*","\\1",colnames(info_dataframe_backup_fictiona_istds))
    ISTD_variabily=t(apply(info_dataframe_backup_fictiona_istds, 1, function(x) tapply(x, colnames(info_dataframe_backup_fictiona_istds), my.median2)))
    
  }
  
  m=as.data.frame(lapply(as.data.frame(ISTD_variabily),my.median2)); ##hacemos media del dataframe de istds para la variabilidad  mas abajo

  
  for(y in rownames(info_dataframe)){
    #y="ND_PLASMA1A"
    #y=rownames(info_dataframe)[1]
    ##calculamos el n de imputaciones realizadas para detectar condiciones a excluir
    #sum(is.na(match(info_dataframe_backup[y,],info_dataframe[y,])))
    #n_imputs=sum(is.na(match(info_dataframe_backup[y,],info_dataframe[y,])))
    n_imputs=which(info_dataframe[y,]!=info_dataframe_backup[y,])
    n_imputs=append(n_imputs,which(is.na(info_dataframe[y,])!=is.na(info_dataframe_backup[y,])))
    n_imputs=unique(n_imputs)
    ##buscamos categorias de las especies imputadas dentro de la variable index
    imputs_category=match(colnames(info_dataframe[,n_imputs]),index[,1],nomatch = 0)
    imputs_category=imputs_category[(imputs_category!=0)]
    imputs_category=table(gsub("^([a-zA-Z]+)(.*)","\\1",index[,2][imputs_category]))
    total_category=table(gsub("^([a-zA-Z]+)(.*)","\\1",index[,2]))
    total_category=total_category[names(imputs_category)]
    
    imputs_category=imputs_category[order(imputs_category,decreasing = T)]
    imputs_category=paste0(names(imputs_category),":",imputs_category,"(",round(imputs_category/total_category[names(imputs_category)]*100),"%);")
    
    n_imputs=length(n_imputs)
    if((n_imputs/ncol(info_dataframe)*100)>porcentaje_resumen){
      resumen_item=paste(if(n_imputs/ncol(info_dataframe)*100>umbral){remove_rows=c(remove_rows,which(rownames(info_dataframe)==rownames(info_dataframe)[which(outlier==T)]));  paste0("->!",format(n_imputs/ncol(info_dataframe)*100,digits=digits),"% de los casos!")} else { " "})
      outliers_resumen_sustituciones[y,"imputaciones"] <-paste("se considera outlier en",n_imputs,"ocasiones de ",ncol(info_dataframe))
      outliers_resumen_sustituciones[y,"comentario"] <-resumen_item
      print(paste("---->",y,"se considera outlier en",n_imputs,"ocasiones de ",ncol(info_dataframe),resumen_item))
      ###buscamos posibles problemas en los ISTD
      rm(t,z)
      t=(as.data.frame(ISTD_variabily[which(rownames(ISTD_variabily)==y),])); ###seleccionamos la fila que presenta problemas
      #a veces me da error porque se genera como un dataframe en lugar de una lista
      if(ncol(t)==1){
        t=t(t)
      }
      z=rbind("outlier"=t,"median"=m);
      z=z[1,]/z[2,]*100-100;
      colnames(z)=gsub("ISTD|IS|_IS|_ISTD","",colnames(z))
      z=z[,colnames(sort(abs(z[1,]),decreasing = T))]
      
      if(n_imputs/ncol(info_dataframe)*100>umbral){
        
        
        ##buscamos todos las muestras de este grupo para valorar posible error durante la normalizacion
        same_group=grep(gsub(paste0("(,*)",regex_filas),"\\1",y),rownames(info_dataframe_backup))
        same_group=same_group[!same_group %in% grep(y,rownames(info_dataframe_backup))] ##buscamos el resto de condiciones del mismo grupo para hacer comparaciones
        same_group=(t(unlist(lapply(info_dataframe_backup[same_group,],my.median))/info_dataframe_backup[y,]))
        
        type_variabilidad=cbind(same_group,gsub("^([a-zA-Z]+)(.*)","\\1",rownames(same_group)))
        type_variabilidad=aggregate(as.numeric(type_variabilidad[,1]),list(type_variabilidad[,2]),my.median)
        colnames(type_variabilidad)=c("Grupo","Factor")
        type_variabilidad=type_variabilidad[order(type_variabilidad[,2],decreasing = T),]
        
        same_group_varibility=same_group[!hboutlier(same_group,r=3)] ### es una medida para quitar un poco de variabilidad al error estandar para el calculo de CV 
        
        if(use_mean==F){
          resumen_sustituciones_median=mean(same_group, na.rm = TRUE)
        } else {
          resumen_sustituciones_median=median(same_group, na.rm = TRUE)
        }
        
        outliers_resumen_sustituciones[y,"factor_cv"] <- (sd(same_group_varibility, na.rm = TRUE)/resumen_sustituciones_median*100)
        outliers_resumen_sustituciones[y,"factor"] <- resumen_sustituciones_median
        if(informe==T){
          print(paste("->Factor de multiplicacion por grupo"))
          print(type_variabilidad)
        }
        print(paste("El factor de multiplicacion medio para la muestra outlier con respecto a su grupo esta en ",format(resumen_sustituciones_median) ))
        print(paste("La dispersion(CV) del factor a lo largo de distintas especies es de ",format((sd(same_group_varibility, na.rm = TRUE)/resumen_sustituciones_median)*100,digits=digits),"%."))
        print(paste("El valor recomendable para aplicarlo es <30%. "))
        
      }
      if(length(z[which(abs(z)>desviacion_ISTD)])>0){
        print(paste0("->ISTDs con desviacion mayor de ",desviacion_ISTD,"% respecto total (la media de clase si no hay ISTDs):"))
        print(paste0( names(z[which(abs(z)>desviacion_ISTD)]),":", round(unlist(z[which(abs(z)>desviacion_ISTD)])),"%; "))
        print(paste0("Clases imputadas, n(%): "))
        print(paste0(imputs_category))
        print(paste0(""))
      }
      print(paste0(""))
    }
  }
  
  n_supera_umral=sum(grepl("!",outliers_resumen_sustituciones$comentario))
  print(paste("---------------------------------------------------------------"))
  print(paste0("--------------->",n_supera_umral," filas superan el umbral establecido en ",umbral,"% <---------------"))
  print(paste("*Se ha creado un backup de los datos originales en la variable info_dataframe_backup"))
  
  if(is.data.frame(index)==T){
    
    # save_imput_images=c("y/n")
    # names(save_imput_images)=c("quieres guardar imagenes de imputacione en formato png? y/n ")
    # fix(save_imput_images)
    # save_imput_images=get("save_imput_images",envir=.GlobalEnv)
    save_imput_images=readline("Quieres guardar imagenes de imputacione en formato png? y/n ")
    if(tolower(str_trim(save_imput_images))=="y"){
      n=0
      dir.create(sprintf(paste0(directory,"/imputations_plots")),showWarnings=F);
      
      while((n*(fig_col*fig_row)-length(unique(as.character(plot_data$NAME))))<0){
        #la grafica propiamente dicha
        n=n+1
        plot_data_new=plot_data[(((fig_col*fig_row)*(n-1)*length(unique(plot_data$ttos)))+1):((fig_col*fig_row)*(n)*length(unique(plot_data$ttos))),]
        
        if(imput_mean==F){
          ### generamos plot distinto si no se hacen imputa
          p<-ggplot(data=plot_data_new, aes(x=grupo, y=value))+ 
            stat_summary(fun.data=function(...) mean_sdl(..., mult=1),geom='errorbar', width=0.5, color='lightblue',size=0.75) + 
            stat_summary(fun.y=my.mean, aes(ymin=..y.., ymax=..y..),geom='errorbar', width=1, color='lightblue', size=1.25)+
            geom_dotplot(binaxis='y', stackdir='center',stackratio=0.01,aes(y=value),dotsize=(dotsize*0.5),fill="black",alpha=0.3)+
            geom_dotplot(position=position_dodge2(0.5), binaxis='y', stackdir='center',stackratio=1,aes(y=value_old),dotsize=dotsize,fill="red",alpha=0.5)+
            geom_label(position=position_jitter(),label.size=0,vjust = 0.2,aes(label=caso,y=value_old, alpha=caso),label.padding = unit(0, "lines"),color="black",size=3)+
            #annotate(geom = "text", x = plot_data[!is.na(plot_data$value_imput),"grupo"], y = plot_data[!is.na(plot_data$value_imput),"value_old"], label = plot_data[!is.na(plot_data$value_imput),"ttos"], size = 4)+
            # lines for means, using geom=errorbar
            theme_bw(15) +
            theme(axis.text.x=element_text(size=8,angle=90,hjust=1.0,vjust=0.5,face="plain"))+
            theme(axis.text.y=element_text(size=8,angle=0,hjust=.5,vjust=.5,face="plain"))
        } else {
          p<-ggplot(data=plot_data_new, aes(x=grupo, y=value))+ 
            stat_summary(fun.data=function(...) mean_sdl(..., mult=1),geom='errorbar', width=0.5, color='lightblue',size=0.75) + 
            stat_summary(fun.y=my.mean, aes(ymin=..y.., ymax=..y..),geom='errorbar', width=1, color='lightblue', size=1.25)+
            geom_dotplot(binaxis='y', stackdir='center',stackratio=0.01,aes(y=value),dotsize=(dotsize*0.5),fill="black",alpha=0.3)+
            geom_dotplot(position=position_dodge2(0.5), binaxis='y', stackdir='center',stackratio=1,aes(y=value_old),dotsize=dotsize,fill="red",alpha=0.5)+
            geom_dotplot(position=position_dodge2(-0.5), binaxis='y', stackdir='center',stackratio=1,aes(y=value_imput),dotsize=dotsize,fill="green",alpha=0.8)+
            geom_label(position=position_jitter(),label.size=0,vjust = 0.2,aes(label=caso,y=value_old, alpha=caso),label.padding = unit(0, "lines"),color="black",size=3)+
            #annotate(geom = "text", x = plot_data[!is.na(plot_data$value_imput),"grupo"], y = plot_data[!is.na(plot_data$value_imput),"value_old"], label = plot_data[!is.na(plot_data$value_imput),"ttos"], size = 4)+
            # lines for means, using geom=errorbar
            theme_bw(15) +
            theme(axis.text.x=element_text(size=8,angle=90,hjust=1.0,vjust=0.5,face="plain"))+
            theme(axis.text.y=element_text(size=8,angle=0,hjust=.5,vjust=.5,face="plain"))
        }
        imput_plot=p + facet_wrap( ~ variable, scales="free", ncol=fig_col, nrow=fig_row)+
          xlab("") +
          #theme(subtitle=)+
          #ylab(paste("Mol% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
          #ylab("nmol/mg")+ #ojo con unidades
          ylab("")+
          #ylab(expression(paste("Mol% of total")))+ #ojo con unidades
          #scale_y_continuous(breaks=0:20*20,expand=c(0,0), limits = c(0, 120))+
          #scale_y_continuous(breaks=0:100*20,expand=c(0,0))+
          #geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
          theme(#legend.position=c(.10, .98), #horizontal, vertical
            legend.position="none",
            strip.background=element_blank(),
            panel.grid.major = element_blank(), # switch off major gridlines
            panel.grid.minor = element_blank(), # switch off minor gridlines
            panel.border = element_blank(),
            axis.line = element_line(colour = "black",size = 0.8, linetype = "solid"),
            axis.ticks.length=unit(0.2, "cm"), 
            #axis.ticks.margin=unit(0.05, "cm"),
            axis.ticks=element_line(size = 0.8, colour = "black", linetype = "solid"),
            #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
            #panel.border = element_border(c("left","bottom")),
            legend.title = element_blank(), # switch off the legend title
            legend.text = element_text(size=12),
            legend.key.size = unit(1.2, "lines")
          ) 
        print(imput_plot)
        
        
        myPng(paste0("imputations_plots/png_imputations_",method,"_r_",Rval,"_n_",n,"_",as.numeric(Sys.time())),16,9, 150);
        print(imput_plot)
        dev.off()
      }
    }
  }
  
  
  ######    APLICAMOS LAS DECISIONES TOMADAS    #########
  

  if(n_supera_umral==0){
    ask_umbral="n"
  } else {

    # ask_umbral=c("y/n/f")
    # names(ask_umbral)=c("Quieres eliminar las condiciones que superan el umbral de % ? y/n (f para factor)")
    # 
    # fix(ask_umbral)
    # ask_umbral=get("ask_umbral",envir=.GlobalEnv)
    ask_umbral=readline("Quieres eliminar las condiciones que superan el umbral de % ? y/n (f para factor)")
    if(tolower(str_trim(ask_umbral))=="f"){
      # umbral_factor=c(35)
      umbral_factor=readline("Indica el valor de CV del factor de multiplicacion maximo(la dispersion maxima) por debajo del cual a las
                             muestra se le aplicara el factor (normalmente 35): ")
      # fix(umbral_factor)
      # umbral_factor=get("umbral_factor",envir=.GlobalEnv)
      outliers_resumen_sustituciones[which(outliers_resumen_sustituciones[,"factor_cv"]>0 & outliers_resumen_sustituciones[,"factor_cv"]>umbral_factor),"factor"]=1
    }
  }
  # ask_imput=c("y/n")
  # names(ask_imput)=c("Quieres aplicar los outliers? y/n ")
  # fix(ask_imput)
  # ask_imput=get("ask_imput",envir=.GlobalEnv)
  ask_imput=readline("Quieres aplicar los outliers? y/n ")
  ##definimos filas para quitar/factorizar
  to_remove=grep("!",outliers_resumen_sustituciones$comentario)
  to_remove=rownames(outliers_resumen_sustituciones)[to_remove]
  
  
  if(tolower(str_trim(ask_imput))=="y"){#si  queremos imputar todos los outliers
    ### introducimos los datos obtenidos en la matriz original
    if(sum(sel_columns!=F)>0 & sum(sel_rows!=F)>0){
      original_info_dataframe[sel_rows,sel_columns]=info_dataframe
      info_dataframe=original_info_dataframe
    } else if(sum(sel_columns!=F)>0 & sum(sel_rows!=F)==0){
      original_info_dataframe[,sel_columns]=info_dataframe
      info_dataframe=original_info_dataframe
    } else if(sum(sel_columns!=F)==0 & sum(sel_rows!=F)>0){
      original_info_dataframe[sel_rows,]=info_dataframe
      info_dataframe=original_info_dataframe
    }
    
    print("!Imputado!")
    to_remove=which(rownames(info_dataframe) %in% to_remove)
    
    ### filtramos las filas si es necesarios
    if(tolower(ask_umbral)=="y" & n_supera_umral!=0){
      info_dataframe=info_dataframe[-to_remove,]
      print("Eliminados los casos outliers!")
    }
    if(tolower(ask_umbral)=="f" & n_supera_umral!=0){
      
      info_dataframe[to_remove,]=info_dataframe_backup[to_remove,]*outliers_resumen_sustituciones[!is.na(outliers_resumen_sustituciones[,"factor"]),"factor"]
      print("Se ha aplicado el factor para los outliers extremos!")
    }
    
    
    sink()
    rm(ask_imput,ask_umbral,umbral_factor)
    return(info_dataframe)
    
  } else { #si no queremos imputar todos los outliers 
    ### introducimos los datos obtenidos en la matriz original
    if(sum(sel_columns!=F)>0 & sum(sel_rows!=F)>0){
      original_info_dataframe[sel_rows,sel_columns]=info_dataframe_backup
      info_dataframe_backup=original_info_dataframe
    } else if(sum(sel_columns!=F)>0 & sum(sel_rows!=F)==0){
      original_info_dataframe[,sel_columns]=info_dataframe_backup
      info_dataframe_backup=original_info_dataframe
    } else if(sum(sel_columns!=F)==0 & sum(sel_rows!=F)>0){
      original_info_dataframe[sel_rows,]=info_dataframe_backup
      info_dataframe_backup=original_info_dataframe
    }
    print("!No imputado!")
    to_remove=which(rownames(info_dataframe_backup) %in% to_remove)
    
    ### filtramos las filas si es necesarios
    if(tolower(ask_umbral)=="y" & n_supera_umral!=0){
      
      info_dataframe_backup=info_dataframe_backup[-to_remove,]
      print("Eliminados los casos outliers!")
    }
    if(tolower(ask_umbral)=="f" & n_supera_umral!=0){
      
      info_dataframe_backup[to_remove,]=info_dataframe_backup[to_remove,]*outliers_resumen_sustituciones[!is.na(outliers_resumen_sustituciones[,"factor"]),"factor"]
      print("Se ha aplicado el factor para los outliers extremos!")
    }
    sink()
    rm(ask_imput,ask_umbral,umbral_factor)
    return((info_dataframe_backup))
  }
  
  
  # # # # FORMATO DE LOS PARAMETROS
  # info_dataframe=t_rsALL.o #matriz de filas para grupos y columnas para especies
  # Rval=2 #2 para hb lo que significa 2 sd o un numero sin porcentaje para CV
  # digits=2 #digitos para los numeros que salen en los mensajes
  # method="hb" #metodo "hb" para doble sd o "CV" para el coeficiente de variacion
  # umbral=30 #umbral para determinar en que porcentaje de especies tiene que imputarse para considerarse una muestra extrema.
  # informe=T #imprimir un informe con cada muestra que ha sido imputada. T o F
  # imput_mean=T #si imputar "T" o eliminar "F" los outliers. influye en las sumas
  # sel_columns=F #si queremos imputas solo algunas columnas ponemos el numero "1:10" para seleccionar una submatriz o FALSE para desactivar
  # sel_rows=F #si queremos imputas solo algunas filas ponemos el numero "1:10" para seleccionar una submatriz o FALSE para desactivar
  # regex_filas="_.*$" #regex para hacer grupos
  # fig_col=4 #columnas para las graficas de outliers imputados
  # fig_row=4 #filas para lo mismo
  # dotsize=3 #tamaÃ±o de los puntos para lo mismo
  # imput_max=F ###numero de valores que quieres imputar porcentaje de muestras que se permiten imputar, enter 0 y 1
  # ISTD_variabily=ISTD_variabilidad ### un dataframe de valores de los ISTD con nombre de filas iguales al los nombres de ttos en info_dataframe o false cuando no lo haya
  # index=rsALL[,c("mz_rt","NAME")] ###coge dos columnas, la primera con nombre de columna de info_dataframe y la segunda su correspondencia para pintar las graficas
  # porcentaje_resumen=5 #para el informe
  # desviacion_ISTD=20 #para el informe
  # NAs_for_median=F #si quieres sustituir todos los NAs por la mediana de su grupo
  # custom_imput_sd=0.01 # es el valor por el cual quieres multiplicar la mediana para disminuir variabilidad durante la imputacion
  # use_mean=F
}
