library(ggpubr)
library(rstatix)

#inputData=t_rsALL.concTOistd.i
index=eindex



#preguntamos numbre de la carpeta para guardar las cosas
directory=get("directory",envir=.GlobalEnv)
print("Estas son las carpetas que hay ahora mismo creadas para analisis:")
print(matrix(gsub("(.*)\\/(.*)$","\\2",list.dirs(directory))))
carpeta_para_figuras=c("figurasPvalue")
names(carpeta_para_figuras) <-c("Introduzca un nombre para la carpeta donde se van a guardar las figuras generadas: ")
fix(carpeta_para_figuras)
#creamos directorio nuevo
dir.create(sprintf(paste0(directory,"/",carpeta_para_figuras)),showWarnings=F);
sub_folder<-paste0(directory,"/",carpeta_para_figuras)
print(paste("Directorio",sub_folder,"creado"))

Basecolor<-c("Blues/Reds/Greens/Yellows/Oranges/Purples")
names(Basecolor) <-c("Elige UNA de las paletas:")
fix(Basecolor)

#########################################
#########################################
# GRAFICOS CON ESTADISTICAS #############
#########################################
#### filtrar ttos
dfc=inputData
dfc$tto=gsub("_\\d{1,}$","",rownames(dfc))
dfc=reshape2::melt(dfc,id=c("tto"))
dfc=summarySE(dfc,measurevar = "value",groupvars = c("variable","tto"))
dfc$variable=gsub("_\\d{1,}$","",dfc$variable)
dfc.sum=ddply(dfc, .(variable,tto), colwise(myfun2,colnames(dfc)[c(4:10)]))
dfc=dfc%>%group_by(variable,tto)%>%summarise(N=mean(N))
dfc.sum=merge(dfc[,1:3],dfc.sum,by=c("variable","tto"))





print("Vamos a hacer un filtrado de datos para las graficas de pvalue")
dfc.sum.sel<-dfc.sum[which(!is.na(dfc.sum$tto)),]

dfc.sum.sel_tto=as.data.frame(cbind(Name=unique(dfc.sum.sel$tto),toRemove=NA))

fix(dfc.sum.sel_tto)
if(length(dfc.sum.sel_tto[!is.na(dfc.sum.sel_tto$toRemove),"Name"])>0){
  dfc.sum.sel<-dfc.sum.sel[(dfc.sum.sel$tto %in% as.character(dfc.sum.sel_tto[is.na(dfc.sum.sel_tto$toRemove),"Name"])),]
}
#### quitar classes

dfc.sum.sel_class=as.data.frame(cbind(class=unique(dfc.sum.sel$variable),toRemove=NA))

fix(dfc.sum.sel_class)
if(length(dfc.sum.sel_class[!is.na(dfc.sum.sel_class$toRemove),"toRemove"])>0){
  dfc.sum.sel<-dfc.sum.sel[(dfc.sum.sel$variable %in% as.character(dfc.sum.sel_class[is.na(dfc.sum.sel_class$toRemove),"class"])),]
}

#dfc.sum.sel$tto<-factor(dfc.sum.sel$tto, levels=unique(dfc.sum.sel$tto)[c(4,1,2,3)], ordered = T)

#### ordenar ttos

dfc.sum.sel_levels=as.data.frame(cbind(levels=unique(dfc.sum.sel$tto),order=NA),stringAsFactor=F)

fix(dfc.sum.sel_levels)
if(length(dfc.sum.sel_levels[!is.na(dfc.sum.sel_levels$order),"order"])>0){
  # which(rownames(dfc.sum.sel_levels) == as.character(dfc.sum.sel_levels$order))
  # 
  # 
  # dfc.sum.sel_levels=dfc.sum.sel_levels[order(dfc.sum.sel_levels$order),]
  # x=as.factor(dfc.sum.sel_levels$levels)
  # y=as.character(dfc.sum.sel_levels$order)
  # factor(x,levels=unique(as.factor(x))[as.numeric(y)],ordered=F)
  dfc.sum.sel$tto <- factor(dfc.sum.sel$tto, levels = as.character(dfc.sum.sel_levels$levels)[match(rownames(dfc.sum.sel_levels),dfc.sum.sel_levels$order)])
  
}

#### seleccionar comparaciones
df <- as.data.frame(matrix(data="",ncol = length(unique(unique(dfc.sum.sel$tto))), nrow = length(unique(unique(dfc.sum.sel$tto)))),stringAsFactor=F)

rownames(df) <- as.character(unique(unique(dfc.sum.sel$tto)))
colnames(df) <- as.character(unique(unique(dfc.sum.sel$tto)))
fix(df)

for(x in colnames(df)){
  #x=colnames(df)[1]
  if(x==colnames(df)[1]){
    comparaciones=list(rows=c(),cols=c())
  }
  rowsToCompare=rownames(df[(as.character(df[,x])!="" & !is.na(df[,x]) & rownames(df)!=x),])
  if(length(rowsToCompare)>0 && rowsToCompare!="NA"){
    comparaciones$rows=c(comparaciones$rows,rowsToCompare)
    comparaciones$cols=c(comparaciones$cols,rep(x,times=length(rowsToCompare)))
  }
}
comparaciones
comparaciones<-mapply(cbind,(comparaciones$rows),(comparaciones$cols),SIMPLIFY = F)

#colnames(dfc.sum.sel)[1]<-"tto"
tabla_sum<-summaryRES_PACK_2(dfc.sum.sel,digitos=3,varTOcast = "tto",estat="media_se",my_comparisons=comparaciones,salto=1.1)
#tabla_sum<-summaryRES_PACK_2(dfc.sum,digitos=3,varTOcast = "tto",estat="media_se",my_comparisons=comparaciones,salto=1.1)
tabla_sum$p_table


dfc.dots=inputData
dfc.dots$tto = gsub("\\..*", "", rownames(dfc.dots))
dfc.dots = reshape2::melt(dfc.dots, id = c("tto"))
dfc.dots$variable = gsub("(\\D+)_(\\d+$)", "\\1", dfc.dots$variable)
dfc.dots = ddply(dfc.dots, .(variable, tto), colwise(myfun2, c("value")))
dfc.dots$tto = gsub("(\\D*)_(\\d+$)", "\\1", dfc.dots$tto)
dfc.sum = merge(dfc.sum.sel,dfc.dots, by = c("variable", "tto"))
colnames(dfc.sum)[which(colnames(dfc.sum) %in% "value.x")]="value"

sub_fusion<-dfc.sum

#sub_fusion<-dfc.sum.sel[which(dfc.sum.sel[,varTOcast] %in% c("WT","12W","30W","CCL4")
#                          & dfc.sum$variable %in% c("Cer","dhCer","SM","dhSM","HexCer","dhHexCer")),]
#& dfc.sum.sel$variable %in% c("CE","FC","TG","PC","PE","AcylCer")),]

##REORDER
#sub_fusion[,varTOcast]<-factor(sub_fusion[,varTOcast],levels=unique(sub_fusion[,varTOcast])[c(2,4,1,3,6,8,5,7)],ordered = T)

###################### TESTING
#escala<-c("darkblue","red","darkred","yellow","darkblue","red","darkred","yellow")
escala=rep(c("white",c(brewer.pal(9,Basecolor)))[c(1,3,5,7,9,10)],times=1)

#escala<-c("darkblue","red","#b35806","#92c5de","#e08214")
#escala<-escala[c(4,1,5,3)]

#puntos<-c(0,15,2,17,5,18)

pd <- position_dodge(0.5) # move them .05 to the left and right

##CLASS BARPLOTS

class_barplot_list<-list();
for(i in unique(as.character(sub_fusion$variable))){
  
  #myPng(paste(paste(i,collapse="_"),"_pvalue_comparison",sep=""),  16,9, 200,dir=sub_folder);
  
  #i="HexCer"
  sub_fusion_for=sub_fusion[which(sub_fusion$variable==i),]
  p_table_sub=tabla_sum$p_table[which(tabla_sum$p_table$supp==i),]
  max_axis_y<-max(p_table_sub$y.position)
  p_table_sub<-p_table_sub[which(!p_table_sub$p.signif %in% c("ns")),]
  
  # The errorbars overlapped, so use position_dodge to move them horizontally
  par(mfcol=c(1,1), mar=c(4, 4, 2, 2), oma=c(0,0,0,0));
  
  
  pd <- position_dodge(1) # move them .05 to the left and right
  # p<-ggplot(sub_fusion_for, aes(x=tto, y=value)) +
  #  geom_bar(stat="identity",aes(fill=tto),color="black")+ 
  # geom_errorbar(aes(ymin=value, ymax=value+sd), width=.5, position=pd,size=0.5) +
  #scale_fill_manual(values=rep (c("white","black"),times=8
  sub_fusion_for$value[-match(unique(sub_fusion_for$value),sub_fusion_for$value)]=NA
  p<-ggplot(sub_fusion_for, aes(x=tto, y=value)) +
    geom_bar(stat="identity",aes(fill=tto,color=tto),size=0.75)+ 
    geom_errorbar(aes(ymin=value, ymax=value+sd,color=tto), width=.5, position=pd,size=0.75,show.legend = F) +
    scale_color_manual(values=rep("black",each=length(sub_fusion_for$tto)))+ #+
  #scale_color_manual(values=rep(escala,each=2))+ #+
  #scale_fill_manual(values=escala)+
    scale_fill_manual(values=escala)
  ## mostrar puntos 
  p=p+geom_dotplot(
    binaxis = 'y',
    stackdir = 'center',
    stackratio = 0.75,
    dotsize = 0.5,
    alpha = 0.8,
    aes(x=tto,y = value.y, fill = tto),
    show.legend = FALSE,
    position = position_dodge(.75)
    
  ) 
  
  auc_plot<-p + #facet_wrap( ~ variable, scales="free",nrow=3,ncol=2)+
    xlab("") +
    ggtitle(i)+
    guides(shape=F,fill=F,color=F)+
    geom_hline(aes(yintercept = 0),
               colour = "black",
               linetype = "solid") +
    #ylab(paste("Mol% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
    ylab("")+ #ojo con unidades
    #ylab(expression(paste("Mol% of total")))+ #ojo con unidades
    theme_bw(13) + 
    theme(plot.title=element_text(size=16,angle=0,hjust=.5,vjust=0.5,face="plain"))+
    theme(axis.text.x=element_text(size=11,angle=0,hjust=.5,vjust=0.5,face="plain"))+
    theme(axis.text.y=element_text(size=12,angle=0,hjust=.5,vjust=.5,face="plain"))+
    scale_y_continuous(expand=c(0,0),limits=c(0,max_axis_y*1.1))+
    scale_x_discrete(labels = c("T2_H_WT" = "CNT","T2_H_12W" = "22w","T2_H_30W" = "30w",
                                "T2_H_CCL4" = "CCl4-22w"))+
    theme(legend.position="top",
          #legend.position=c(.9, .10),
          #legend.position = "none",
          strip.background=element_blank(),
          panel.grid.major = element_blank(), # switch off major gridlines
          panel.grid.minor = element_blank(), # switch off minor gridlines
          panel.border = element_blank(),
          axis.line = element_line(),
         # axis.text.x=element_text(size=12),
          #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
          #panel.border = element_border(c("left","bottom")),
          legend.title = element_blank(), # switch off the legend title
          legend.text = element_text(size=8),
          legend.key.size = unit(1.0, "lines") 
          
    ) 
  if(all(p_table_sub$p.signif %in% c("ns"))==F){
    auc_plot<- auc_plot + stat_pvalue_manual(
      p_table_sub,xmin="group1",xmax = "group2", size=3,
      label = "p.signif",tip.length = 0.01,
      #position = position_dodge(0.8), ## hace que las brackets se vean mal
      bracket.size = 0.5,
      bracket.nudge.y = 0,
      #bracket.shorten=0.5,
      step.increase=0.02,## para a?adir espacio entre las brackets
      hide.ns=T
    )
  }
  class_barplot_list[[i]]<-auc_plot
  # ggsave(filename = paste0(directory,"/graficos/AUC_ALL_",i,".png"), auc_plot_list[[i]],
  #        width = 7.5, height = 8, dpi = 600, units = "cm", device='png')
  
  #print(auc_plot)
  
  #pushViewport(viewport(y=0.98,x=0.85,height=.05))
  #title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
  #grid.table(title)
  
  #dev.off()
  
}

##CLASS BOXPLOTS

class_boxplot_list<-list();
for(i in unique(as.character(sub_fusion$variable))){
  
  #myPng(paste(paste(i,collapse="_"),"_pvalue_dotplot_comparison",sep=""),  16,9, 200,dir=sub_folder);
  #i="SM"
  sub_fusion_for=sub_fusion[which(sub_fusion$variable==i),]
  p_table_sub=tabla_sum$p_table[which(tabla_sum$p_table$supp==i),]
  p_table_sub<-p_table_sub[which(!p_table_sub$p.signif %in% c("ns")),]
  
  # The errorbars overlapped, so use position_dodge to move them horizontally
  par(mfcol=c(1,1), mar=c(4, 4, 2, 2), oma=c(0,0,0,0));
  
  
  pd <- position_dodge(1) # move them .05 to the left and right
  # p<-ggplot(sub_fusion_for, aes(x=tto, y=value)) +
  #  geom_bar(stat="identity",aes(fill=tto),color="black")+ 
  # geom_errorbar(aes(ymin=value, ymax=value+sd), width=.5, position=pd,size=0.5) +
  #scale_fill_manual(values=rep (c("white","black"),times=8
  sub_fusion_for$value[-match(unique(sub_fusion_for$value),sub_fusion_for$value)]=NA
  
  p<- 
    ggplot(sub_fusion_for, aes(x=tto, y=value.y)) +
    stat_boxplot(geom = "errorbar", width = 0.4,size=0.75)+
    geom_boxplot(outlier.shape = NA, show.legend = F,alpha=1,notch=F,lwd=0.5,aes(x=tto, y=value.y, fill=tto))+
    geom_beeswarm(aes(x=tto, y=value.y, color=tto),
                  cex=1.5,
                  stroke=.5,
                  size=.5,
                  alpha=.75, 
                  position=position_dodge(0.0), 
                  binaxis='y', 
                  stackdir='center',
                  stackratio=1, 
                  show.legend = F)+
    # geom_dotplot(aes(x=tto, y=value.y, fill=tto, color=tto),binaxis = 'y',stackdir = 'center',stackratio = 0.75,dotsize = 0.5,alpha = 0.8, 
    #              show.legend = FALSE,position = position_dodge(.75))+
    # 
    scale_fill_manual(values=escala)+
    scale_color_manual(values=escala[order(escala,decreasing=F)])
  
  auc_plot<-p + #facet_wrap( ~ variable, scales="free",nrow=3,ncol=2)+
    xlab("") +
    ggtitle(i)+
    guides(shape=F,fill=F,color=F)+
    #geom_hline(aes(yintercept = 0),
    #          colour = "black",
    #         linetype = "solid") +
    #ylab(paste("Mol% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
    ylab("")+ #ojo con unidades
    #ylab(expression(paste("Mol% of total")))+ #ojo con unidades
    theme_bw(13) + 
    theme(plot.title=element_text(size=16,angle=0,hjust=.5,vjust=0.5,face="plain"))+
    theme(axis.text.x=element_text(size=11,angle=0,hjust=.5,vjust=0.5,face="plain"))+
    theme(axis.text.y=element_text(size=12,angle=0,hjust=.5,vjust=.5,face="plain"))+
    #scale_y_continuous(expand=c(0,0),limits=c(0,max(p_table_sub$y.position)*1.2))+
    #scale_x_discrete(labels = c("T2_H_WT" = "CNT","T2_H_12W" = "22w","T2_H_30W" = "30w",
    #"T2_H_CCL4" = "CCl4-22w"))+
    theme(legend.position="top",
          #legend.position=c(.9, .10),
          #legend.position = "none",
          strip.background=element_blank(),
          panel.grid.major = element_blank(), # switch off major gridlines
          panel.grid.minor = element_blank(), # switch off minor gridlines
          panel.border = element_blank(),
          axis.line = element_line(),
          # axis.text.x=element_text(size=12),
          #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
          #panel.border = element_border(c("left","bottom")),
          legend.title = element_blank(), # switch off the legend title
          legend.text = element_text(size=8),
          legend.key.size = unit(1.0, "lines") 
          
    ) 
  if(all(p_table_sub$p.signif %in% c("ns"))==F){
    auc_plot<- auc_plot + stat_pvalue_manual(
      p_table_sub,xmin="group1",xmax = "group2", size=3,
      label = "p.signif",tip.length = 0.01,
      #position = position_dodge(0.8), ## hace que las brackets se vean mal
      bracket.size = 0.5,
      bracket.nudge.y = 0,
      #bracket.shorten=0.5,
      step.increase=0.02, ## para a?adir espacio entre las brackets
      hide.ns=T
    )
  }
  class_boxplot_list[[i]]<-auc_plot
  # ggsave(filename = paste0(directory,"/graficos/AUC_ALL_",i,".png"), auc_plot_list[[i]],
  #        width = 7.5, height = 8, dpi = 600, units = "cm", device='png')
  
  #print(auc_plot)
  
  #pushViewport(viewport(y=0.98,x=0.85,height=.05))
  #title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
  #grid.table(title)
  
  #dev.off()
  
}


####GAFICOS DE COMPARACIONES ENTRE ESPECIES####

puntitos=c("y/n")
names(puntitos) <-c("Quieres pintar puntos: ")
fix(puntitos)

dfc=inputData
dfc$tto=gsub("(.*)_(\\d+$)","\\1",rownames(dfc))
dfc=reshape2::melt(dfc,id=c("tto"))
dfc=as.data.frame(summarySE(dfc,measurevar = "value",groupvars = c("tto","variable")))
dfc$TAG=gsub("(^.*)([\\s_ ])([^\\s_]+)$","\\1",dfc$variable)
dfc$m_NAME=rep(index[(index$index %in% dfc$variable),"NAME"],times=length(unique(dfc$tto)))

tabla_esp<-summaryRES_PACK_2(dfc,digitos=3,varTOcast = "tto",estat="media_se",my_comparisons=comparaciones,salto=1.1)
tabla_esp$stat<-cbind(NAME=index[match(tabla_esp$stat$variable,index$index),"NAME"],tabla_esp$stat)
# cambio de nombres
# dfc$tto<-revalue(dfc$tto,c(revalue_plotlegent))
# dfc$tto<-factor(dfc$tto,levels=plotlegend_custom)

# ### mostrar puntos
dfc.dots=inputData
dfc.dots$tto = gsub("\\..*", "", rownames(dfc.dots))
dfc.dots = reshape2::melt(dfc.dots, id = c("tto"))
dfc.dots = ddply(dfc.dots, .(variable, tto), colwise(myfun2, c("value")))
dfc.dots$tto = gsub("(.*)_(\\d+$)", "\\1", dfc.dots$tto)
dfc = merge(dfc,dfc.dots, by = c("variable", "tto"))
colnames(dfc)[which(colnames(dfc) %in% "value.x")]="value"
# 
unique(dfc$tto)
###transform to long format
dfc.long<-dfc
head(dfc.long, 3)

dfc.long<-dfc.long[(dfc.long$tto %in% as.character(dfc.sum.sel_tto[is.na(dfc.sum.sel_tto$toRemove),"Name"])),]
dfc.long<-dfc.long[(dfc.long$TAG %in% as.character(dfc.sum.sel_class[is.na(dfc.sum.sel_class$toRemove),"class"])),]
dfc.long$tto <- factor(dfc.long$tto, levels = as.character(dfc.sum.sel_levels$levels)[match(rownames(dfc.sum.sel_levels),dfc.sum.sel_levels$order)])

###Correct FC class just to respresent the specie with CE
dfc.long[which(dfc.long$TAG %in% "FC"),"TAG"]<-"CE"


#custom_order<-c("NORMAL","NAFL","NASH","NASHFSIG")
escala=rep(c("white",c(brewer.pal(9,Basecolor)))[c(1,3,5,7,9,10)],times=3)

#dfc.long$tto<-factor(dfc.long$tto,levels=custom_order,ordered=T)
referencia<-c(as.character(unique(dfc.long$tto)))
referencia=selectOption(referencia,name="select_reference_sample")

step.upper<-1.15 # salto de eje-y para que queden bien las barras de comparaciones
#!is.na(referencia)
species_barplot_list<-list()
#stat_test_barplot_list<-list()
region<-unique(dfc.long$TAG)
for (cls in region){
  
  #cls<-"HexCer"
  #creamos un subset con regiones indicadas arriba
  sub_fusion<-subset(dfc.long, dfc.long$TAG %in% cls)
  #saltamos grupos sin datos
  if(nrow(sub_fusion)>0){
    #quitamos los plasmalogenos de las graficas
    sub_fusion=sub_fusion[!grepl("O-",sub_fusion$m_NAME),]
    
    #seleccionamos las 10 especies mas abundantes
    sub_fusion_new=c()
    for( x in unique(sub_fusion$TAG)){
      top_ten=names(sort(tapply(sub_fusion[which(sub_fusion$TAG==x),"value"], sub_fusion[which(sub_fusion$TAG==x),"variable"], max),decreasing = T))[1:9]
      sub_fusion_new=rbind(sub_fusion_new,sub_fusion[which(sub_fusion$variable %in% top_ten),])
    }
    sub_fusion=sub_fusion_new
    
    sub_fusion$TAG<-factor(sub_fusion$TAG)
    
    
    ## quitamos nombres de grupo/region excluyendo FC para que no se quede sin nombre region
    if(!(cls %in% c("FC","DhS1P","S1P","dhSphC18","SphC18","HexSph"))){
      sub_fusion$m_NAME<-as.character(gsub(paste0(paste0(cls,"| |\\*+")), "",sub_fusion$m_NAME)) 
    }
    sub_fusion$variable=index[match(sub_fusion$variable,index$index),"NAME"]
    
    par(mfcol=c(1,1))
    
    sub_fusion$tto
    ##añado pq falla en algunos casos el ajuste del barplot tengo que ver pq.
    sub_fusion$value.y[sub_fusion$value.y==0]<-NA
    ###
    sub_fusion$m_NAME<-factor(sub_fusion$m_NAME,levels=sort(unique(sub_fusion$m_NAME)),ordered=T)
    
    #statistical test
    stat.test <- sub_fusion %>%
      group_by(m_NAME) %>%
      t_test(value.y ~ tto) %>%
      adjust_pvalue(method = "none") %>%
      add_significance("p.adj")
    #stat.test<-stat.test[which(!stat.test$p.adj.signif %in% c("ns")),]
    #x=factor(m_NAME,levels=sort(unique(m_NAME)),ordered=T) 
    # Create a bar plot with error bars (mean +/- sd)
    bp <- ggbarplot(
      sub_fusion, x = "m_NAME", y = "value.y", 
      #add = "mean_sd", error.plot = "upper_errorbar",
      add="mean",
      fill= "tto", palette = escala , 
      width=0.8, #anchura de barra
      size=0.5, #tamaño borde de barras
      position = position_dodge(.8)
    )+geom_errorbar(aes(ymin=value,ymax=value+sd,fill=tto),
                    width=.5,
                    size=.5,
                    position=position_dodge(.8))
    
    # Add p-values onto the bar plots
    if(length(referencia)!=0){
      stat.test <- stat.test %>%
        add_xy_position(fun = "mean_sd", x = "m_NAME", dodge = 0.8, ref.group=referencia,step.increase = 0.05) 
    }else{
      stat.test <- stat.test %>%
        add_xy_position(fun = "mean_sd", x = "m_NAME", dodge = 0.8,step.increase = 0.05,comparisons=comparaciones) 
    }
    # bp + stat_pvalue_manual(
    #   stat.test,  label = "p.adj.signif", tip.length = 0.01
    # )
    # 
    # Move down the brackets using `bracket.nudge.y`
    bp + stat_pvalue_manual(
      stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns=T,
      bracket.nudge.y = -2
    )
    
    bp<-bp + stat_pvalue_manual(
      stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns=T
      #remove.bracket = T, vjust = 2
    )
    if(tolower(puntitos)=="y"){
        bp<-bp+geom_dotplot(
          binaxis = 'y',
          stackdir = 'center',
          stackratio = 0.8,
          dotsize = .4,
          alpha = 0.8,
          aes(x=m_NAME,y = value.y, fill = tto),
          #aes(x=m_NAME,y = value.y),
          show.legend = FALSE,
          position = position_dodge(.8)
        )
    }
    bp<-bp + #facet_wrap( ~ variable, scales="free",ncol=3)+ ##Coment to separate classes
      xlab("") +
      #ylab(paste("% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
      ylab("") +
      #geom_hline(aes(yintercept = 0),
      #          colour = "black",
      #         linetype = "solid") +
      theme_classic(base_line_size = 1) +
      labs(title = cls) +
      scale_y_continuous(expand=c(0,0),limits=c(0,max(stat.test$y.position)*step.upper))+
      #scale_y_continuous(expand=c(0,0))+
      #scale_x_continuous(expand = expansion(mult = c(0, .2))) +
      #guides( fill = FALSE)+#labs(subtitle = x)+
      #guides(shape=F)+
      theme(
        axis.line = element_line(colour = 'black', size = 0.5), ##control line width of x and y axis
        axis.ticks = element_line(colour = 'black', size = 0.5), ##control thick width of x and y tips
        #axis.text.x = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = rel(1.2)),
        axis.text.y = element_text(
          size = 12,
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
        #legend.position="top",
        #legend.position=c(.9, .10),
        legend.position = "none",
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
    species_barplot_list[[cls]]<-bp
    #stat_test_barplot_list[[cls]]<-stat.test
    print("?Barplot de especies hecho!")
    
    #print(bp)
    
    #pushViewport(viewport(y=0.95,x=0.85,height=.05))
    #title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
    #grid.table(title)
    
    #dev.off()
  }
  
}

# #species_plot_list
# dfc.long.sel<-dfc.long[,c(1:10)]
# #colnames(dfc.long.sel)[1]<-"variable"
# dfc.long.sel$variable<-as.character(dfc.long.sel$variable)
# dfc.long.sel<-dfc.long.sel[!duplicated(dfc.long.sel$variable),]
# tabla_esp<-summaryRES_PACK_2(dfc.long.sel[c(1:3),],digitos=4,varTOcast = "tto",estat="media_se",my_comparisons=comparaciones,salto=1.1)


#### BOXPLOT SPECIES
species_boxplot_list<-list()
stat_test_boxplot_list<-list()
for (cls in region){
  
  #cls<-"TG"
  #creamos un subset con regiones indicadas arriba
  sub_fusion<-subset(dfc.long, dfc.long$TAG %in% cls)
  #saltamos grupos sin datos
  if(nrow(sub_fusion)>0){
    #quitamos los plasmalogenos de las graficas
    sub_fusion=sub_fusion[!grepl("O-",sub_fusion$m_NAME),]
    
    #seleccionamos las 10 especies mas abundantes
    sub_fusion_new=c()
    for( x in unique(sub_fusion$TAG)){
      top_ten=names(sort(tapply(sub_fusion[which(sub_fusion$TAG==x),"value"], sub_fusion[which(sub_fusion$TAG==x),"variable"], max),decreasing = T))[1:9]
      sub_fusion_new=rbind(sub_fusion_new,sub_fusion[which(sub_fusion$variable %in% top_ten),])
    }
    sub_fusion=sub_fusion_new
    
    sub_fusion$TAG<-factor(sub_fusion$TAG)
    
    
    ## quitamos nombres de grupo/region excluyendo FC para que no se quede sin nombre region
    if(!(cls %in% c("FC","DhS1P","S1P","dhSphC18","SphC18","HexSph"))){
      sub_fusion$m_NAME<-as.character(gsub(paste0(paste0(cls,"| |\\*+")), "",sub_fusion$m_NAME)) 
    }
    sub_fusion$variable=index[match(sub_fusion$variable,index$index),"NAME"]
    
    par(mfcol=c(1,1))
    
    sub_fusion$tto
    ##añado pq falla en algunos casos el ajuste del barplot tengo que ver pq.
    sub_fusion$value.y[sub_fusion$value.y==0]<-NA
    ###
    sub_fusion$m_NAME<-factor(sub_fusion$m_NAME,levels=sort(unique(sub_fusion$m_NAME)),ordered=T)
    
    #statistical test
    stat.test <- sub_fusion %>%
      group_by(m_NAME) %>%
      t_test(value.y ~ tto) %>%
      adjust_pvalue(method = "none") %>%
      add_significance("p.adj")
    
    bp<-ggboxplot(
      sub_fusion, x = "m_NAME", y = "value.y", 
      fill= "tto", 
      palette = escala, 
      #width=0.8, #anchura de barra
      #size=0.5,
      outlier.shape = NA, #tamaño borde de barras
      lwd=0.25
      #position = position_dodge(.8)
    ) 
    
    # Add p-values onto the bar plots
    if(length(referencia)!=0){
      stat.test <- stat.test %>%
        add_xy_position(fun = "median_iqr", x = "m_NAME", dodge = 0.8, ref.group=referencia,step.increase = 0.1) 
    }else{
      stat.test <- stat.test %>%
        add_xy_position(fun = "median_iqr", x = "m_NAME", dodge = 0.8, step.increase = 0.1,comparisons=comparaciones) 
    }
    bp + stat_pvalue_manual(
      stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns=T, label.size =0.25 ,
      bracket.nudge.y = -2
    )
    
    bp<-bp + stat_pvalue_manual(
      stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns=T
      #remove.bracket = T, vjust = 2
    )
    if(tolower(puntitos)=="y"){
      bp<-bp+geom_dotplot(
        binaxis = 'y',
        stackdir = 'center',
        stackratio = 0.1,
        dotsize = .1,
        alpha = 0.8,
        aes(x=m_NAME,y = value.y, fill = tto),
        #aes(x=m_NAME,y = value.y),
        show.legend = FALSE,
        position = position_dodge(.8)
      )
    }
    bp<-bp + #facet_wrap( ~ variable, scales="free",ncol=3)+ ##Coment to separate classes
      xlab("") +
      #ylab(paste("% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
      ylab("") +
      #geom_hline(aes(yintercept = 0),
      #          colour = "black",
      #         linetype = "solid") +
      theme_classic(base_line_size = 1) +
      labs(title = cls) +
      scale_y_continuous(expand=c(0,0),limits=c(0,max(stat.test$y.position)*step.upper))+
      #scale_y_continuous(expand=c(0,0))+
      #scale_x_continuous(expand = expansion(mult = c(0, .2))) +
      #guides( fill = FALSE)+#labs(subtitle = x)+
      #guides(shape=F)+
      theme(
        axis.line = element_line(colour = 'black', size = 0.5), ##control line width of x and y axis
        axis.ticks = element_line(colour = 'black', size = 0.5), ##control thick width of x and y tips
        #axis.text.x = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = rel(1.2)),
        axis.text.y = element_text(
          size = 12,
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
        #legend.position=c(.9, .10),
        legend.position = "top",
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
    
    species_boxplot_list[[cls]]<-bp
    #stat_test_boxplot_list[[cls]]<-stat.test
    print("Boxplot de especies hecho!!!")
    
    #print(bp)
    
    #pushViewport(viewport(y=0.95,x=0.85,height=.05))
    #title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
    #grid.table(title)
    
    #dev.off()
  }
  
}

###GRAFICOS CLASES
ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_class_Boxplot_sphingos.png"),
       annotate_figure(
         ggarrange(class_boxplot_list$dhCer, class_boxplot_list$Cer,
                   class_boxplot_list$dhSM, class_boxplot_list$SM,
                   class_boxplot_list$dhHexCer, class_boxplot_list$HexCer,
                   #labels = c("A", "B", "C","D","E","F"),
                   ncol = 2, nrow = 3),
         top = text_grob(paste0("Quantitation of lipid classes in ", select_type), color = "black", face = "plain", size = 14,just="center")
       ),
       width = 20, height = 22, dpi = 600, units = "cm", device='png')

ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_class_Boxplot_neutral.png"),
       annotate_figure(
         ggarrange(class_boxplot_list$CE, class_boxplot_list$FC,
                   class_boxplot_list$TG, class_boxplot_list$ACer,
                   class_boxplot_list$PE, class_boxplot_list$LPC,
                   #labels = c("A", "B", "C","D","E","F"),
                   ncol = 2, nrow = 3),
         top = text_grob(paste0("Quantitation of lipid classes in ", select_type), color = "black", face = "plain", size = 14,just="center")
       ),
       width = 20, height = 22, dpi = 600, units = "cm", device='png')

# ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_class_Barplot_sphingos.png"),
#        annotate_figure(
#          ggarrange(class_barplot_list$dhCer, class_barplot_list$Cer,
#                    class_barplot_list$dhSM, class_barplot_list$SM,
#                    class_barplot_list$dhHexCer, class_barplot_list$HexCer,
#                    #labels = c("A", "B", "C","D","E","F"),
#                    ncol = 2, nrow = 3),
#          top = text_grob(paste0("Quantitation of lipid classes in ", select_type), color = "black", face = "plain", size = 14,just="center")
#        ),
#        width = 20, height = 22, dpi = 600, units = "cm", device='png')
# 
# ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_class_Barplot_neutral.png"),
#        annotate_figure(
#          ggarrange(class_barplot_list$CE, class_barplot_list$FC,
#                    class_barplot_list$TG, class_barplot_list$ACer,
#                    class_barplot_list$PC, class_barplot_list$PE,
#                    #labels = c("A", "B", "C","D","E","F"),
#                    ncol = 2, nrow = 3),
#          top = text_grob(paste0("Quantitation of lipid classes in ", select_type), color = "black", face = "plain", size = 14,just="center")
#          ),
#        width = 20, height = 22, dpi = 600, units = "cm", device='png')


###GRAFICOS ESPECIES
ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_species_Boxplot_neutral.png"),
       annotate_figure(
         ggarrange(species_boxplot_list$CE, species_boxplot_list$TG,
                   species_boxplot_list$ACer, species_boxplot_list$PE,
                   species_boxplot_list$PC, species_boxplot_list$LPC,
                   #labels = c("A", "B", "C","D","E","F"),
                   common.legend = TRUE,
                   ncol = 2, nrow = 3),
         top = text_grob(paste0("Quantitation of lipid species in ", select_type), color = "black", face = "plain", size = 14,just="center")
       ),
       width = 20, height = 22, dpi = 600, units = "cm", device='png')

ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_species_Boxplot_sphingos.png"),
       annotate_figure(
         ggarrange(species_boxplot_list$dhCer, species_boxplot_list$Cer,
                   species_boxplot_list$dhSM, species_boxplot_list$SM,
                   species_boxplot_list$dhHexCer, species_boxplot_list$HexCer,
                   #labels = c("A", "B", "C","D","E","F"),
                   common.legend = TRUE,
                   ncol = 2, nrow = 3),
         #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
         top = text_grob(paste0("Quantitation of lipid species in ", select_type), color = "black", face = "plain", size = 14,just="center")
       ),
       width = 20, height = 22, dpi = 600, units = "cm", device='png')

