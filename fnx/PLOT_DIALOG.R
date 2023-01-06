  #inputData=t_rsALL.concTOistd.i
  # eindex=rsALL.samples[,c("eindex","NAME_ltr")]
  eindex<-data.frame(cbind("index"=colnames(rsALL),"NAME"=colnames(rsALL)))
  
  #preguntamos numbre de la carpeta para guardar las cosas
  directory=get("directory",envir=.GlobalEnv)
  print("Estas son las carpetas que hay ahora mismo creadas para analisis:")
  print(matrix(gsub("(.*)\\/(.*)$","\\3",list.dirs(directory))))
  carpeta_para_figuras=c("figuras")
  names(carpeta_para_figuras) <-c("Introduzca un nombre para la carpeta donde se van a guardar las figuras generadas: ")
  fix(carpeta_para_figuras)
  #creamos directorio nuevo
  dir.create(sprintf(paste0(directory,"/",carpeta_para_figuras)),showWarnings=F);
  sub_folder<-paste0(directory,"/",carpeta_para_figuras)
  print(paste("Directorio",sub_folder,"creado"))

  inputData_tto=unique(gsub("(.*)_(\\d+$)","\\1",rownames(inputData)))
  inputData_tto=as.data.frame(cbind(Name=unique(inputData_tto),toRemove=NA))
  
  fix(inputData_tto)
  if(length(inputData_tto[!is.na(inputData_tto$toRemove),"Name"])>0){
    
    #paste0("^",as.character(inputData_tto[is.na(inputData_tto$toRemove),"Name"]),collapse="|")
    inputData<-inputData[grepl(paste0("^",as.character(inputData_tto[is.na(inputData_tto$toRemove),"Name"]),collapse="|"),rownames(inputData)),]
  }
  #### quitar classes
  inputData_class=unique(gsub("(.*)_(\\d+$)","\\1",colnames(inputData)))
  inputData_class=as.data.frame(cbind(class=unique(inputData_class),toRemove=NA))
  
  fix(inputData_class)
  if(length(inputData_class[!is.na(inputData_class$toRemove),"toRemove"])>0){
    inputData<-inputData[grepl(paste0("^",as.character(inputData_class[is.na(inputData_class$toRemove),"class"]),collapse = "|"),colnames(inputData))]
  }
  
 

  ####---- Haciendo graficas nmol/mg---------###
  
   
    
    
    dfc=inputData
    dfc$tto=gsub("(.*)_(\\d+$)","\\1",rownames(dfc))
    dfc=reshape2::melt(dfc,id=c("tto"))
    dfc=summarySE(dfc,measurevar = "value",groupvars = c("variable","tto"))
    dfc$variable=gsub("_\\d{1,}$","",dfc$variable)
    dfc=ddply(dfc, .(variable,tto), colwise(myfun2,colnames(dfc)[c(4:10)]))
    
    
    # cambio de nombres
    #dfc$tto<-revalue(dfc$tto,c(revalue_plotlegend))
    # 
    # dfc$tto<-factor(dfc$tto,levels=plotlegend_custom)
    ## variables generales
    
    plotlegend_custom=unique(dfc[,"tto"], incomparables = FALSE, fromLast = TRUE) 
    pallete=colorRampPalette(c("white","grey"))(length(plotlegend_custom))
    region=unique(gsub("_\\d{1,}","",colnames(inputData)))
    
    
    
    dfc.dots=inputData
    dfc.dots$tto = gsub("\\..*", "", rownames(dfc.dots))
    dfc.dots = reshape2::melt(dfc.dots, id = c("tto"))
    dfc.dots$variable = gsub("(.*)_(\\d+$)", "\\1", dfc.dots$variable)
    dfc.dots = ddply(dfc.dots, .(variable, tto), colwise(myfun2, c("value")))
    dfc.dots$ID=dfc.dots$tto
    dfc.dots$tto = gsub("(.*)_(\\d+$)", "\\1", dfc.dots$tto)
    dfc = merge(dfc,dfc.dots, by = c("variable", "tto"))
    colnames(dfc)[which(colnames(dfc) %in% "value.x")]="value"
    
    ## resultados en tabla para exportar por si acaso
    class_barplot_data=dfc
    
    #### Barplot de clases
    for(cls in region) {
      #cls=region[1]
      sub_fusion<-subset(dfc, dfc$variable %in% cls & dfc$tto %in% plotlegend_custom)
      
      #saltamos grupos sin datos
      if(nrow(sub_fusion)>0){
        
        
        sort(unique(sub_fusion$variable))
        p <-ggplot(data = sub_fusion, aes(x=variable,
                                          y=value,
                                          fill=factor(tto,levels=plotlegend_custom,ordered=T)))+
          geom_bar(position=position_dodge(),stat="identity",width=0.75,size=0.75) +
          geom_bar(position=position_dodge(),stat="identity",color="black",width=.75,size=.75,show.legend=FALSE) + #to get rid of stripped legend 
          
          geom_errorbar(aes(ymin=value-sd,ymax=value+sd),
                        width=.4,
                        size=.75,
                        position=position_dodge(.75))+
          scale_fill_manual(name = "", values = as.character(pallete), #change colours
                            labels =plotlegend_custom
          )+
          scale_colour_continuous(guide = F)
        ## mostrar puntos 
        p=p+geom_dotplot(
          binaxis = 'y',
          stackdir = 'center',
          stackratio = 0.8,
          dotsize = .6,
          alpha = 0.8,
          aes(x=variable,y = value.y, fill = tto),
          show.legend = FALSE,
          position = position_dodge(.75)
          
        ) 
        
        
        
        myPng(paste("sum_barplot_",paste(cls,collapse="_"),sep=""),  16,9, 200,dir=sub_folder);
        
        
        p= p +# facet_wrap( ~ TAG, scales="free",nrow=2,ncol=2)+
          xlab("") +
          #ylab(paste("% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
          ylab(custom_units_to_use)+ #ojo con unidades
          #ylab("Area/mg Prote?na")+ #ojo con unidades
          geom_hline(aes(yintercept = 0),
                     colour = "black",
                     linetype = "solid") +
          theme_classic(base_line_size = 1) +
          labs(title = cls) +
          scale_y_continuous(expand = expansion(mult = c(0, .2))) +
          #guides( fill = FALSE)+#labs(subtitle = x)+
          #guides(shape=F)+
          theme(
            axis.text.x = element_blank(),
            plot.caption = element_text(hjust = 0.5, size = rel(1.2)),
            axis.text.y = element_text(
              size = 16,
              angle = 0,
              hjust = .5,
              vjust = .5,
              face = "plain",
              color = "black"
            ),
            plot.title = element_text(hjust = 0.5, size = 16),
            axis.title.y = element_text(size = 16, margin = margin(
              t = 0,
              r = 20,
              b = 0,
              l = 0
            )),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            #legend.position="right",
            #legend.position=c(.9, .10),
            #legend.position = "none",
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            # switch off major gridlines
            panel.grid.minor = element_blank(),
            # switch off minor gridlines
            panel.border = element_blank(),
            #axis.line = element_line(),
            #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
            #panel.border = element_border(c("left","bottom")),
            legend.title = element_blank(),
            # switch off the legend title
            legend.text = element_text(size = 16),
            legend.key.size = unit(1.0, "lines")
            
          )
        print(p)
        
        pushViewport(viewport(y=0.95,x=0.85,height=.05))
        title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
        grid.table(title)
        
        dev.off()
      }
    }
    print("?Barplot de clases hecho!")
    
    ##Barplot de especies
    # 
    # dfc=inputData
    # dfc$tto=gsub("(.*)_(\\d+$)","\\1",rownames(dfc))
    # dfc=reshape2::melt(dfc,id=c("tto"))
    # dfc=as.data.frame(summarySE(dfc,measurevar = "value",groupvars = c("tto","variable")))
    # dfc$TAG=gsub("(^.*)([\\s_ ])([^\\s_]+)$","\\1",dfc$variable)
    # dfc$m_NAME=rep(eindex[(eindex$index %in% dfc$variable),"NAME"],times=length(unique(dfc$tto)))
    # 
    # 
    # 
    # # cambio de nombres
    # # dfc$tto<-revalue(dfc$tto,c(revalue_plotlegent))
    # # dfc$tto<-factor(dfc$tto,levels=plotlegend_custom)
    # 
    # # ### mostrar puntos
    # dfc.dots=inputData
    # dfc.dots$tto = gsub("\\..*", "", rownames(dfc.dots))
    # dfc.dots = reshape2::melt(dfc.dots, id = c("tto"))
    # dfc.dots = ddply(dfc.dots, .(variable, tto), colwise(myfun2, c("value")))
    # dfc.dots$tto = gsub("(.*)_(\\d+$)", "\\1", dfc.dots$tto)
    # dfc = merge(dfc,dfc.dots, by = c("variable", "tto"))
    # colnames(dfc)[which(colnames(dfc) %in% "value.x")]="value"
    # # 
    # species_barplot_data=dfc
    # 
    # for (cls in region) {
    #   #creamos un subset con regiones indicadas arriba
    #   sub_fusion<-subset(dfc, dfc$TAG %in% cls & dfc$tto %in% plotlegend_custom)
    #   
    #   #saltamos grupos sin datos
    #   if(nrow(sub_fusion)>0){
    #     #quitamos los plasmalogenos de las graficas
    #     sub_fusion=sub_fusion[!grepl("O-",sub_fusion$m_NAME),]
    #     
    #     #seleccionamos las 10 especies mas abundantes
    #     sub_fusion_new=c()
    #     for( x in unique(sub_fusion$TAG)){
    #       top_ten=names(sort(tapply(sub_fusion[which(sub_fusion$TAG==x),"value"], sub_fusion[which(sub_fusion$TAG==x),"variable"], max),decreasing = T))[1:9]
    #       sub_fusion_new=rbind(sub_fusion_new,sub_fusion[which(sub_fusion$variable %in% top_ten),])
    #     }
    #     sub_fusion=sub_fusion_new
    #     
    #     sub_fusion$TAG<-factor(sub_fusion$TAG)
    #     
    #     
    #     ## quitamos nombres de grupo/region excluyendo FC para que no se quede sin nombre region
    #     if(!(cls %in% c("FC","DhS1P","S1P","dhSphC18","SphC18","HexSph"))){
    #       sub_fusion$m_NAME<-as.character(gsub(paste0(paste0(cls,"| |\\*+")), "",sub_fusion$m_NAME)) 
    #     }
    #     sub_fusion$variable=eindex[match(sub_fusion$variable,eindex$index),"NAME"]
    #     
    #     par(mfcol=c(1,1))
    #     
    #     
    #     p <-ggplot(data = sub_fusion, aes(x=factor(m_NAME,levels=sort(unique(m_NAME)),ordered=T),
    #                                       y=value,
    #                                       fill=factor(tto,levels=plotlegend_custom,ordered=T))) + 
    #       geom_bar(position=position_dodge(),stat="identity",width=0.75,size=0.75) +
    #       geom_bar(position=position_dodge(),stat="identity",color="black",width=.75,size=.75,show.legend=FALSE) + #to get rid of stripped legend 
    #       
    #       geom_errorbar(aes(ymin=value,ymax=value+sd),
    #                     width=.4,
    #                     size=.75,
    #                     position=position_dodge(.75))+
    #       scale_fill_manual(name = "", values = as.character(pallete), #change colours
    #                         labels = plotlegend_custom)+
    #       scale_colour_continuous(guide = F) #+
    #     
    #     # ## mostrar puntos 
    #      p=p+geom_dotplot(
    #        binaxis = 'y',
    #        stackdir = 'center',
    #       stackratio = 0.8,
    #       dotsize = .6,
    #       alpha = 0.8,
    #       aes(x=m_NAME,y = value.y, fill = tto),
    #       show.legend = FALSE,
    #       position = position_dodge(.75)
    #     )
    #     
    #     myPng(paste("molecule_barplot_",paste(cls,collapse="_"),sep=""),  16,9, 200,dir=sub_folder);
    #     
    #     p = p + facet_wrap( ~ variable, scales="free",ncol=3)+
    #       xlab("") +
    #       #ylab(paste("% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
    #       ylab(custom_units_to_use) +
    #       geom_hline(aes(yintercept = 0),
    #                  colour = "black",
    #                  linetype = "solid") +
    #       theme_classic(base_line_size = 1) +
    #       labs(title = cls) +
    #       scale_y_continuous(expand = expansion(mult = c(0, .2))) +
    #       #guides( fill = FALSE)+#labs(subtitle = x)+
    #       #guides(shape=F)+
    #       theme(
    #         axis.text.x = element_blank(),
    #         plot.caption = element_text(hjust = 0.5, size = rel(1.2)),
    #         axis.text.y = element_text(
    #           size = 16,
    #           angle = 0,
    #           hjust = .5,
    #           vjust = .5,
    #           face = "plain",
    #           color = "black"
    #         ),
    #         plot.title = element_text(hjust = 0.5, size = 16),
    #         axis.title.y = element_text(size = 16, margin = margin(
    #           t = 0,
    #           r = 20,
    #           b = 0,
    #           l = 0
    #         )),
    #         axis.title.x = element_blank(),
    #         axis.ticks.x = element_blank(),
    #         #legend.position="right",
    #         #legend.position=c(.9, .10),
    #         #legend.position = "none",
    #         strip.background = element_blank(),
    #         panel.grid.major = element_blank(),
    #         # switch off major gridlines
    #         panel.grid.minor = element_blank(),
    #         # switch off minor gridlines
    #         panel.border = element_blank(),
    #         #axis.line = element_line(),
    #         #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
    #         #panel.border = element_border(c("left","bottom")),
    #         legend.title = element_blank(),
    #         # switch off the legend title
    #         legend.text = element_text(size = 16),
    #         legend.key.size = unit(1.0, "lines")
    #         
    #       ) 
    #     print(p)
    #     
    #     pushViewport(viewport(y=0.95,x=0.85,height=.05))
    #     title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
    #     grid.table(title)
    #     
    #     dev.off()
    #   }
    #   
    # }
    # print("?Barplot de especies hecho!")
    ##### dotplots for all molecules
    
    # dfc=inputData
    # dfc$tto=rownames(dfc)
    # dfc$tto=gsub("_\\d{1,}$","",rownames(dfc))
    # dfc=reshape2::melt(dfc,id=c("tto"))
    # dfc=summarySE(dfc,measurevar = "value",groupvars = c("tto","variable"))
    # dfc$TAG=gsub("(^.*)([\\s_ ])([^\\s_]+)$","\\1",dfc$variable)
    # dfc$m_NAME=rep(eindex[(eindex$index %in% dfc$variable),"NAME"],times=length(unique(dfc$tto)))
    # dfc$m_cond=gsub("(.*)\\D(\\d+$)","\\2\\3",dfc$tto)
    # dfc$tto=gsub("(.*)\\D(\\d+$)","\\1",dfc$tto)
    # 
    
    dfc=inputData
    dfc$tto=rownames(dfc)
    dfc=reshape2::melt(dfc,id=c("tto"))
    dfc=summarySE(dfc,measurevar = "value",groupvars = c("tto","variable"))
    dfc$TAG=gsub("(^.*)([\\s_ ])([^\\s_]+)$","\\1",dfc$variable)
    dfc$m_NAME=rep(eindex[(eindex$index %in% dfc$variable),"NAME"],times=length(unique(dfc$tto)))
    dfc$m_cond=gsub("(.*)_(\\d+$)","\\2\\3",dfc$tto)
    dfc$tto=gsub("(.*)_(\\d+$)","\\1",dfc$tto)
    
    
    
    
    # cambio de nombres 
    # dfc$tto<-revalue(dfc$tto,c(revalue_plotlegent))
    # 
    # dfc$tto<-factor(dfc$tto,levels=plotlegend_custom)
    
    species_dotplot_data=dfc
    for (cls in region) {
      #cls<-"dhCer"
      #creamos un subset con regiones indicadas arriba
      #sub_fusion<-subset(dfc.sum, dfc.sum$variable %in% cls & dfc.sum$tto %in% plotlegend)
      sub_fusion<-subset(dfc, dfc$TAG %in% cls & dfc$tto %in% plotlegend_custom)
      
      #saltamos grupos sin datos
      if(nrow(sub_fusion)>0){
        
        # #quitamos los plasmalogenos de las graficas
        # sub_fusion=sub_fusion[!grepl("O-",sub_fusion$m_NAME),]
        
        #seleccionamos las 10 especies mas abundantes
        sub_fusion_new=c()
        for( x in unique(sub_fusion$TAG)){
          top_ten=names(sort(tapply(sub_fusion[which(sub_fusion$TAG==x),"value"], sub_fusion[which(sub_fusion$TAG==x),"variable"], my.max),decreasing = T))[1:9]
          sub_fusion_new=rbind(sub_fusion_new,sub_fusion[which(sub_fusion$variable %in% top_ten),])
        }
        sub_fusion=sub_fusion_new
        
        
        
        #quitamos los sufijos en los nombres
        sub_fusion$m_NAME=gsub("#.*","",sub_fusion$m_NAME)
        
        
        sub_fusion$TAG<-factor(sub_fusion$TAG)
        
        
        sub_fusion$TAG<-factor(sub_fusion$TAG,levels=cls,ordered=T)
        
        sub_fusion$tto=factor(sub_fusion$tto,levels=(plotlegend_custom),ordered=T)
        sub_fusion=sub_fusion[order(sub_fusion$tto),]
        
        ## quitamos nombres de grupo/region excluyendo FC para que no se quede sin nombre region
        if(!(cls %in% c("FC","DhS1P","S1P","dhSphC18","SphC18","HexSph"))){
          sub_fusion$m_NAME<-as.character(gsub(paste0(paste0(cls,"| |\\*+")), "",sub_fusion$m_NAME)) 
        }
        
        
        par(mfcol=c(1,1))
        
        myPng(paste("molecule_dotplot_",paste(cls,collapse="_"),sep=""),  16,9, 200,dir=sub_folder);
        
        p <- ggplot(sub_fusion, aes(x=tto, y=value, fill=(factor(tto,levels=plotlegend_custom,ordered=T)))) + 
          scale_fill_manual(labels = plotlegend_custom, name ="", values=rep(pallete,length(unique(sub_fusion$m_NAME))))
        p = p + theme(legend.position = "top") + 
          geom_boxplot(aes(y=value), width=1) +
          geom_dotplot(binaxis='y', stackdir='center',stackratio=0.01,dotsize=0.5,aes(y=value, fill=tto)) +
          geom_label(label.size=0,position=position_dodge2(0.5),aes(label=m_cond,y=value, alpha=m_cond),label.padding = unit(0, "lines"),color="black",size=3)+
          facet_wrap( ~ TAG*m_NAME, scales="free",ncol=3)+ #nrow=ceil(i/2),ncol=2)+
          xlab("") + 
          #ylab(paste("% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
          ylab(custom_units_to_use)+ #ojo con unidades
          #ylab("Area/mg Prote?na")+ #ojo con unidades
          
          #ylab(expression(paste("% of total")))+ #ojo con unidades
          theme_bw(15) +
          theme(axis.text.x=element_text(size=8,angle=90,hjust=1.0,vjust=0.5,face="plain"))+
          theme(axis.text.y=element_text(size=8,angle=0,hjust=.5,vjust=.5,face="plain"))+
          scale_y_continuous(expand = expansion(mult = c(0, .2))) +
          #scale_y_continuous(breaks=0:10*0.2,expand=c(0,0),limits = c(0, 1.1))+
          geom_hline(aes(yintercept=0), colour="black", linetype="solid")+
          theme(#legend.position=c(.63, .9),
            legend.position="none",
            strip.background=element_blank(),
            panel.grid.major = element_blank(), # switch off major gridlines
            panel.grid.minor = element_blank(), # switch off minor gridlines
            panel.border = element_blank(),
            axis.line = element_line(),
            #axis.text.x=element_text(size=12,angle=45,hjust=1.0,vjust=1.,face="plain"),
            #axis.text.x=element_text(size=12,angle=0,hjust=1.,vjust=1.0,face="plain"),
            #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
            #panel.border = element_border(c("left","bottom")),
            legend.title = element_blank(), # switch off the legend title
            legend.text = element_text(size=12),
            legend.key.size = unit(1.0, "lines") 
          ) 
        
        
        print(p)
        
        pushViewport(viewport(y=0.98,x=0.85,height=.05))
        title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
        grid.table(title)
        
        dev.off()
      }
      
    }
    print("?Dotplot de especies hecho!")
    
  
    
    
    dfc=inputData
    
    dfc$tto=rownames(dfc)
    dfc=reshape2::melt(dfc,id=c("tto"))
    dfc=summarySE(dfc,measurevar = "value",groupvars = c("variable","tto"))
    dfc$variable=gsub("(^.*)([\\s_ ])([^\\s_]+)$","\\1",dfc$variable)
    dfc=ddply(dfc, .(variable,tto), colwise(myfun2,colnames(dfc)[c(4:10)]))
    dfc$TAG=dfc$variable
    dfc$m_cond=gsub("(.*)_(\\d+$)","\\2\\3",dfc$tto)
    dfc$tto=gsub("(.*)_(\\d+$)","\\1",dfc$tto)

    # cambio de nombres
    # dfc$tto<-revalue(dfc$tto,c(revalue_plotlegend))
    # dfc$tto<-factor(dfc$tto,levels=plotlegend_custom)
    
    class_dotplot_data=dfc
    
    for (cls in region) {
      
      sub_fusion<-subset(dfc, dfc$TAG %in% cls & dfc$tto %in% plotlegend_custom)  
      
      #saltamos grupos sin datos
      if(nrow(sub_fusion)>0){
        
        
        
        sub_fusion$TAG<-factor(sub_fusion$TAG)
        
        
        sub_fusion$TAG<-factor(sub_fusion$TAG,levels=cls,ordered=T)
        
        sub_fusion$tto=factor(sub_fusion$tto,levels=(plotlegend_custom),ordered=T)
        sub_fusion=sub_fusion[order(sub_fusion$tto),]
        
        
        
        par(mfcol=c(1,1))
        
        
        myPng(paste("sum_dotplot_",paste(cls,collapse="_"),sep=""),  16,9, 200,dir=sub_folder);
        
        p <- ggplot(sub_fusion, aes(x=tto, y=value, fill=(factor(tto,levels=plotlegend_custom,ordered=T)))) + 
          scale_fill_manual(labels = plotlegend_custom, name ="", values=rep(pallete,length(unique(sub_fusion$TAG))))
        
        p = p + theme(legend.position = "top") + 
          geom_boxplot(aes(y=value), width=1) +
          geom_dotplot(binaxis='y', stackdir='center',stackratio=0.01,dotsize=0.5,aes(y=value, fill=tto)) +
          geom_label(label.size=0,position=position_dodge2(0.5),aes(label=m_cond,y=value, alpha=m_cond),label.padding = unit(0, "lines"),color="black",size=6)+
          #facet_wrap( ~ TAG*m_NAME, scales="free",ncol=4,nrow=3)+ #nrow=ceil(i/2),ncol=2)+
          xlab("") + 
          #ylab(paste("% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
          ylab(custom_units_to_use)+ #ojo con unidades
          #ylab("Area/mg Prote?na")+ #ojo con unidades
          #ylab(expression(paste("% of total")))+ #ojo con unidades
          theme_bw(15) +
          theme(axis.text.x=element_text(size=18,angle=90,hjust=1.0,vjust=0.5,face="plain"))+
          theme(axis.text.y=element_text(size=18,angle=0,hjust=.5,vjust=.5,face="plain"))+
          labs(caption=unique(sub_fusion$TAG)) + 
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
          scale_y_continuous(expand = expansion(mult = c(0, .2))) +
          #scale_y_continuous(breaks=0:10*0.2,expand=c(0,0),limits = c(0, 1.1))+
          geom_hline(aes(yintercept=0), colour="black", linetype="solid")+
          theme(#legend.position=c(.63, .9),
            legend.position="none",
            strip.background=element_blank(),
            panel.grid.major = element_blank(), # switch off major gridlines
            panel.grid.minor = element_blank(), # switch off minor gridlines
            panel.border = element_blank(),
            axis.line = element_line(),
            #axis.text.x=element_text(size=12,angle=45,hjust=1.0,vjust=1.,face="plain"),
            #axis.text.x=element_text(size=12,angle=0,hjust=1.,vjust=1.0,face="plain"),
            #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
            #panel.border = element_border(c("left","bottom")),
            legend.title = element_blank(), # switch off the legend title
            legend.text = element_text(size=12),
            legend.key.size = unit(1.0, "lines") 
          ) 
        
        
        print(p)
        
        pushViewport(viewport(y=0.98,x=0.85,height=.05))
        title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
        grid.table(title)
        
        dev.off()
      }
      
    }
    print("?Dotplot de clases hecho!")
    #save(list = ls(all.names = TRUE)[grepl("backup",ls(all.names = TRUE))], file = "plot_parameters.rda", envir = .GlobalEnv)
    
  ### guardamos los valores usados para generar las graficas
     get_plot_legend_results=list("class_barplot_data"=class_barplot_data,
                 "species_barplot_data"=species_barplot_data,
                 "class_dotplot_data"=class_dotplot_data,
                 "species_dotplot_data"=species_dotplot_data)
     
     #View(get_plot_legend_results)
     
    

