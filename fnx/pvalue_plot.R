library(ggpubr)
library(rstatix)

#inputData=t_rsALL.concTOistd.i
#index = cbind("index"=colnames(rsALL),"NAME"=colnames(rsALL))



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



#########################################
#########################################
# GRAFICOS CON ESTADISTICAS #############
#########################################
#### filtrar ttos
dfc=inputData
dfc$tto=str_split_fixed(rownames(dfc),"#",2)[,1]
dfc=reshape2::melt(dfc,id=c("tto"))
dfc=summarySE(dfc,measurevar = "value",groupvars = c("variable","tto"))
#dfc$variable=gsub("_\\d{1,}$","",dfc$variable)
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

dfc.sum.sel_class=as.data.frame(cbind(class=as.character(unique(dfc.sum.sel$variable)),toRemove=NA))

fix(dfc.sum.sel_class)
if(length(dfc.sum.sel_class[!is.na(dfc.sum.sel_class$toRemove),"toRemove"])>0){
  dfc.sum.sel<-dfc.sum.sel[(dfc.sum.sel$variable %in% as.character(dfc.sum.sel_class[is.na(dfc.sum.sel_class$toRemove),"class"])),]
}


#### ordenar ttos

dfc.sum.sel_levels=as.data.frame(cbind(levels=unique(dfc.sum.sel$tto),order=NA),stringAsFactor=F)

fix(dfc.sum.sel_levels)
if(length(dfc.sum.sel_levels[!is.na(dfc.sum.sel_levels$order),"order"])>0){
  
  dfc.sum.sel$tto <- factor(dfc.sum.sel$tto, levels = as.character(dfc.sum.sel_levels$levels)[match(rownames(dfc.sum.sel_levels),dfc.sum.sel_levels$order)])
  
}

#### seleccionar comparaciones
df <- as.data.frame(matrix(data="",ncol = length(unique(unique(dfc.sum.sel$tto))), nrow = length(unique(unique(dfc.sum.sel$tto)))),stringAsFactor=F)

rownames(df) <- as.character(rev(unique(dfc.sum.sel$tto)))
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
tabla_sum<-summaryRES_PACK_2(dfc.sum.sel,digitos=3,varTOcast = "tto",estat="media_sd",my_comparisons=comparaciones,salto=1)
tabla_sum$p_table


dfc.dots=inputData
dfc.dots$tto = rownames(dfc.dots)
dfc.dots = reshape2::melt(dfc.dots, id = c("tto"))
#dfc.dots$variable = gsub("(\\D+)_(\\d+$)", "\\1", dfc.dots$variable)
#dfc.dots = ddply(dfc.dots, .(variable, tto), colwise(myfun2, c("value")))
#dfc.dots = ddply(dfc.dots, .(variable, tto), colwise(function(x)sum(x[(x > 0) & !is.na(x)],na.rm=T), c("value")))


dfc.dots$tto = str_split_fixed(dfc.dots$tto,"#",2)[,1]
dfc.sum = merge(dfc.sum.sel,dfc.dots, by = c("variable", "tto"))
colnames(dfc.sum)[which(colnames(dfc.sum) %in% "value.x")]="value"

sub_fusion<-dfc.sum


###################### TESTING

n_groups=length(unique(str_split_fixed(sub_fusion$tto,"->",2)[,1]))
n_hk=length(unique(str_split_fixed(sub_fusion$tto,"->",2)[,2]))

escala=rep(c(brewer.pal(9,"Set3"))[1:n_groups],n_hk)



pd <- position_dodge(0.5) # move them .05 to the left and right

library(ggpubr)

## formula aproximada inventada para ajustar el tama?o de etiquetas en el eje x, 
#se divide el 100 entre n de grupos y luego se reduce un porcentaje segun la longitud maxima de las etiquetas
#x_axis_font_size=round((100/length(unique(str_split_fixed(sub_fusion$tto,"#",2)[,1])))*(1-(max(nchar(as.character(sub_fusion$tto)))/100)))


x_axis_font_size<-6
auc_plot_list<-list();
for(i in unique(as.character(sub_fusion$variable))){
 
  #myPng(paste(paste(i,collapse="_"),"_pvalue_comparison",sep=""),  16,9, 200,dir=sub_folder);
  
  #i<-"DEGS1"
  sub_fusion_for=sub_fusion[which(sub_fusion$variable==i),]
  p_table_sub=tabla_sum$p_table[which(tabla_sum$p_table$supp==i),]
  
  
  # The errorbars overlapped, so use position_dodge to move them horizontally
  par(mfcol=c(1,1), mar=c(4, 4, 2, 2), oma=c(0,0,0,0));
  
  
  pd <- position_dodge(1) # move them .05 to the left and right
  # p<-ggplot(sub_fusion_for, aes(x=tto, y=value)) +
  #  geom_bar(stat="identity",aes(fill=tto),color="black")+ 
  # geom_errorbar(aes(ymin=value, ymax=value+sd), width=.5, position=pd,size=0.5) +
  #scale_fill_manual(values=rep (c("white","black"),times=8
  match(unique(sub_fusion_for$tto),sub_fusion_for$tto)
  sub_fusion_for$value[-unique(match(unique(sub_fusion_for$tto),sub_fusion_for$tto),match(unique(sub_fusion_for$value),sub_fusion_for$value))]=NA
  #provisional para evitar puntos que no hay
  #sub_fusion_for[which(sub_fusion_for$value.y==0),"value.y"]<-NA
  
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
    stackratio = 0.8,
    dotsize = 1,
    alpha = 0.8,
    aes(x=tto,y = value.y, fill = tto),
    show.legend = F,
    #na.rm=T
    position = position_dodge(1)
    
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
    theme(axis.text.x=element_text(size=x_axis_font_size,angle=0,hjust=.5,vjust=0.5,face="plain"))+
    theme(axis.text.y=element_text(size=12,angle=0,hjust=.5,vjust=.5,face="plain"))+
    scale_y_continuous(expand=c(0,0),limits=c(0,max(p_table_sub$y.position)*1.5))+
    #scale_x_discrete(labels = c("T2_H_WT" = "CNT","T2_H_12W" = "22w","T2_H_30W" = "30w",
    #                            "T2_H_CCL4" = "CCl4-22w"))+
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
          
    ) + stat_pvalue_manual(
      p_table_sub,xmin="group1",xmax = "group2", size=8,
      label = "p.signif",tip.length = 0.01,
      #position = position_dodge(0.8), ## hace que las brackets se vean mal
      bracket.size = 0.5,
      bracket.nudge.y = 0,
      #bracket.shorten=0.5,
      step.increase=0.02 ## para a?adir espacio entre las brackets
      #hide.ns=T
    )
  
  auc_plot_list[[i]]<-auc_plot
  # ggsave(filename = paste0(directory,"/graficos/AUC_ALL_",i,".png"), auc_plot_list[[i]],
  #        width = 7.5, height = 8, dpi = 600, units = "cm", device='png')
  
  print(auc_plot)
  
   pushViewport(viewport(y=0.98,x=0.85,height=.05))
   title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
   grid.table(title)

   dev.off()
  
}
auc_plot_list
###GRAFICO AGRUPADO
#dir.create(sprintf(paste0(directory,"/graficosTest")),showWarnings=F);
# ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/","qPCR.png"),
#        annotate_figure(
#          ggarrange(plotlist = auc_plot_list,
#                    #labels = c("A", "B", "C","D","E","F"),
#                    ncol = 3, nrow = 3)#,
#          #top = text_grob("Quantitation relative to LTR and IS (not in LTR)", color = "black", face = "bold", size = 16)
#        ),
#        width = 15, height = 22, dpi = 600, units = "cm", device='png')
ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_qPCR.png"),
       annotate_figure(
         ggarrange(plotlist = auc_plot_list,
                   ncol = 4, nrow = 4,
                   common.legend = TRUE),
         top = text_grob(paste0("qPCR in ", select_type), color = "black", face = "plain", size = 20,just="center")
       ),
       width = 30, height = 30, dpi = 300, units = "cm", device='png')




