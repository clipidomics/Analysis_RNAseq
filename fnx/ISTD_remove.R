

#vemos los estandares y la concentracion con sus respectivas variaciones
ISTDs_todos=rsALL[grep(paste0(rsISTD$NAME,collapse="|"),str_trim(rsALL$NAME)),][,c(10,which(colnames(rsALL)==rsSample$Muestra[1]):length(rsALL))]
rownames(ISTDs_todos)=ISTDs_todos$NAME
ISTDs_todos=ISTDs_todos[,-1]

rsALL[which(str_trim(rsALL$NAME) %in% str_trim(rsISTD$NAME)),"NAME"]

 
#mostramos la gr?fica
ISTDs_todos$class=rownames(ISTDs_todos)
ISTD_variabilidad=reshape2::melt(ISTDs_todos)
ISTD_variabilidad$variable=sub("(.*)_\\d+$","\\1",ISTD_variabilidad$variable)
ISTD_variabilidad=summarySE(ISTD_variabilidad,"value",groupvars = "class")[,c(1:3,7)]
ISTD_variabilidad=cbind(ISTD_variabilidad,"Incluir"="y")



### Comprobamos si hay clases con mas de un ISTD  
ISTD_match=(table(rsALL[which(rsALL$type=="ISTD"),"Region"]))

if(length(ISTD_match[which(ISTD_match!=1)])>0){

  if(sum(grepl(0,ISTD_match))==0){
      ##hacemos una figura para comprobar si hay clases con varios ISTDs
    istd_plot=reshape2::melt(ISTDs_todos)
    istd_means=aggregate(istd_plot,by=list(istd_plot$class),FUN = mean)[,c(1,4)]
    colnames(istd_means)=c("class","mean")
    istd_plot=merge(istd_plot,istd_means,id="class")
    colnames(istd_plot)[2]=colnames(rsSample)[1]
    istd_plot=merge(istd_plot,rsSample[,1:2],id="muestra")
    istd_plot$class=sub(" 18:1_17:0","",istd_plot$class)
    
    # creamos una carpeta
    dir.create(sprintf(paste0(getwd(),"/figuras_ISTDs")),showWarnings=F);
    #hacemos la gráfica
    for(cls in names(ISTD_match[which(ISTD_match!=1)])){
      istd_plot_sub=istd_plot[grepl(cls,istd_plot$class),]
      
      p<-ggplot(istd_plot_sub, aes(x=Muestra, y=(value), fill=class)) +
        geom_bar(position=position_dodge(0.5),width=0.75,size=0.75, stat = "identity")+
        geom_bar(aes(y=mean),width=0.75,size=0.75, stat = "identity",fill="black",alpha = 0.2)+
        scale_shape_manual(values=c(0,15,1,16,2,17,5,18))
      #
      istd_plot_p=p + facet_wrap( ~ class+Tanda, scales="free_x",ncol=length(unique(rsSample$Tanda)))+
        xlab("Gray bar represents global mean for each ISTD") +
        ylab("Intensity")+ #ojo con unidades
        theme_bw(13) +
        theme(axis.text.x=element_text(size=5,angle=90,hjust=.5,vjust=.5,face="plain"))+
        theme(axis.text.y=element_text(size=8,angle=90,hjust=.5,vjust=.5,face="plain"))+
        #scale_x_continous(breaks = c(0,4,8,12,16), labels = c("0","4","8","12","16"))+
        #scale_x_continuous(breaks=c(0,4,8,12,16),expand=c(0,0), limits = c(-.5, 17))+
        #scale_y_continuous(expand=c(1,2))+
        theme(
          legend.position = "none",
          strip.background=element_blank(),
          panel.grid.major = element_blank(),  #switch off major gridlines
          panel.grid.minor = element_blank(),  #switch off minor gridlines
          panel.border = element_blank(),
          axis.line = element_line(),
          
        )
      print(istd_plot_p)
      ##guardamos la imagen
      print(paste("Se ha generado un grafico de la calidad de los ISTDS duplicados para:",cls))
      myPng(paste0("/figuras_ISTDs/ISTD_variabilidad_barplot",cls), 16,9, 200);
      print(istd_plot_p)
      dev.off()
    }
  } else {
    
    stop(paste("Falta ISTD de",
               paste0(names(ISTD_match[grepl(0,ISTD_match)]),collapse=",")
               ))
  }
  
  
}

if (sum(ISTD_match)==length(ISTD_match)) {
  print("OK - Hay estandares para cada especie") 
} else { 
  ISTD_not_match = ISTD_match[which(ISTD_match!=median(ISTD_match))]
  print(paste("Problema con siguiente ISTD:",names(ISTD_not_match),"n=",ISTD_not_match, collapse=" "))
  print("Elige qué ISTDs quieremos eliminar del análisis");
  fix(ISTD_variabilidad);
  if(sum(ISTD_variabilidad$Incluir=="n")>0){
    rsALL=rsALL[-grep(paste(ISTD_variabilidad[which(ISTD_variabilidad$Incluir=="n"),"class"],collapse="|"),rsALL$NAME),]
  }
  
} 


### quitar istd que no estan
rsISTD=rsISTD[which(rsISTD$NAME %in% rsALL$NAME),]
if(sum((rsALL$type=="ISTD" & !(as.character(rsALL$Region) %in% as.character(unique(rsALL[rsALL$type!="ISTD","Region"])))))>0){
  rsALL=rsALL[-which(rsALL$type=="ISTD" & !(as.character(rsALL$Region) %in% as.character(unique(rsALL[rsALL$type!="ISTD","Region"])))),]
  rsALL$Region=factor(rsALL$Region)
}
rm(list=ls(all.names = T)[grepl("^istd",ls(all.names = T),ignore.case = T)])
