
input_temp=reshape2::dcast(input, Muestra~Gen,value.var = "Mean.Cp")

#### ANALISIS QC-control para del istd
t_rsALL.istd_rsd<-data.frame(input_temp[,c("Muestra",unlist(unique(input[which(input$Type=="Housekeeping"),"Gen"])))])
rownames(t_rsALL.istd_rsd)=t_rsALL.istd_rsd$Muestra


### TODO OBTENER compatibilizar con analisis lipidomico
GRUPO<-input[match(t_rsALL.istd_rsd$Muestra,input$Muestra),"GRUPO"]
SAMPLEid<-input[match(t_rsALL.istd_rsd$Muestra,input$Muestra),"Muestra"]
colnames(SAMPLEid)="SAMPLEid"
tanda<-input[match(t_rsALL.istd_rsd$Muestra,input$Muestra),"Analysis"]
colnames(tanda)="tanda"

## quitamos la columna de muestras por compatibilidad
t_rsALL.istd_rsd=t_rsALL.istd_rsd[-1]


### continuamos como siempre
t_rsALL.istd_rsd<-cbind(GRUPO,SAMPLEid,tanda,t_rsALL.istd_rsd)

t_rsALL.istd_rsd$GRUPO<-as.character(t_rsALL.istd_rsd$GRUPO)
# t_rsALL.istd_rsd[which(t_rsALL.istd_rsd$GRUPO %in% c(select_calname)),"GRUPO"]<-"QC"
# t_rsALL.istd_rsd[which(!t_rsALL.istd_rsd$GRUPO %in% c(select_calname)),"GRUPO"]<-"SAMPLE"
t_rsALL.istd_rsd$GRUPO<-paste(t_rsALL.istd_rsd$GRUPO,t_rsALL.istd_rsd$tanda,sep="#")
row_ids<-rownames(t_rsALL.istd_rsd)



rownames(t_rsALL.istd_rsd)<-paste(t_rsALL.istd_rsd$GRUPO,t_rsALL.istd_rsd$SAMPLEid,sep="_")
t_rsALL.istd_rsd.i=t_rsALL.istd_rsd[,-c(1:3)]

# t_rsALL.istd_rsd.i<-new_outlier(t_rsALL.istd_rsd[,-c(1:3)],Rval=10,method="hb",imput_mean = T,umbral = 30,informe = T, 
#                                 regex_filas="_.*$",fig_col=4,fig_row=3,dotsize = 1.5,index=t_rsALL.istd_rsd[,c(1:2)],
#                                 ISTD_variabily=F, NAs_for_median=T,
#                                 custom_imput_sd=0.01,
#                                 use_mean=T)

t_rsALL.istd_rsd.i<-cbind(t_rsALL.istd_rsd[,c(1:3)],t_rsALL.istd_rsd.i)

t_rsALL.istd_rsd.i<-t_rsALL.istd_rsd.i[order(t_rsALL.istd_rsd.i$tanda,t_rsALL.istd_rsd.i$GRUPO),]

t_rsALL.istd_rsd.i<-cbind(Run_Number=1:nrow(t_rsALL.istd_rsd.i),t_rsALL.istd_rsd.i)

#match(t_rsALL.istd_rsd.i$SAMPLEid,rsSample$Muestra)
#reodenar según orden de inyección
#TODO

### CONSTRUIR GRAFICO VARIABILIDAD ISTDs

library(ggQC)

m_rsALL.istd_rsd.i<-reshape2::melt(t_rsALL.istd_rsd.i,id=c("GRUPO","tanda","SAMPLEid","Run_Number"),measure=colnames(t_rsALL.istd_rsd.i)[-c(1:4)]);
m_rsALL.istd_rsd.i$variable<-gsub("(.*)\\_(.*)","\\1",m_rsALL.istd_rsd.i$variable)
m_rsALL.istd_rsd.i$GRUPO<-factor(m_rsALL.istd_rsd.i$GRUPO, levels=unique(m_rsALL.istd_rsd.i$GRUPO), ordered=T)

### corregimos la escala
m_rsALL.istd_rsd.i=(m_rsALL.istd_rsd.i%>%group_by(variable)%>%summarise(GRUPO=ifelse(grepl("QC",GRUPO),"QC","Sample"),tanda=tanda,SAMPLEid=SAMPLEid,Run_Number=Run_Number,variable=variable,value=value,valueScale = (value/(10^(nchar(gsub("\\..*","",my.median(value)))-2)))))


#sel.istd<-m_rsALL.istd_rsd.i[which(m_rsALL.istd_rsd.i$variable %in% c("SM")),]
#escala<-c("#d73027","black","#1a9850","black","#3288bd","black")
escala<-rep(c("#d73027","black"),times=100)
dir.create(sprintf(paste0(directory,"/QC_shewhart_plot")),showWarnings=F);
sel_ISTD=list()
for(x in unique(m_rsALL.istd_rsd.i$variable)){
  sel.istd<-m_rsALL.istd_rsd.i[which(m_rsALL.istd_rsd.i$variable %in% x),]
  show_names=1
  sel_ISTD[[x]]=sel.istd
  XmR_Plot <-
    ggplot(sel.istd, aes(x = Run_Number, y = valueScale)) +
    geom_point((aes(fill = GRUPO)),
               size = 3,
               shape = 21,
               colour = "black",
               stroke = 0.5
    ) + geom_line() + # add the points and lines
    stat_QC(
      method = "XmR",
      # specify QC charting method
      auto.label = T,
      # Use Autolabels
      label.digits = 0,
      # Use two digit in the label
      show.1n2.sigma = F   # Show 1 and two sigma lines
    ) +
    
    scale_x_continuous(expand =  expand_scale(mult = .25)) +
    scale_y_continuous(expand =  expand_scale(mult = .15)) +
    scale_fill_manual(values = escala)
  
  ##### dispersion dots
  if(show_names==1){
    XmR_Plot = XmR_Plot + geom_label(
      label.size = 0,
      position = position_dodge2(0.5),
      aes(
        #label = gsub("^([^\\s_]+).*([\\s_]).*([^\\s_]+)$", "\\3", SAMPLEid, perl = T),
        label = Run_Number,
        y = valueScale + (valueScale * 0.15)
      ),
      label.padding = unit(0.1, "lines"),
      label.r	= unit(0.5, "lines"),
      color = "black",
      size = 4,
      show.legend = FALSE
    )
  }
  XmR_Plot2 <-
    XmR_Plot + facet_wrap(~ variable * tanda,
                          scales = "free_x",
                          ncol = length(unique(sel.istd$tanda))) +
    xlab("") +
    ylab(paste0("Area/Intensity x ",sel.istd$value[1]/sel.istd$valueScale[1])) + #ojo con unidades
    #ylab(expression(paste("% of total")))+ #ojo con unidades
    theme_bw(15) +
    theme(
      #axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      #legend.position=c(.63, .9),
      legend.position = "none",
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      # switch off major gridlines
      panel.grid.minor = element_blank(),
      # switch off minor gridlines
      panel.border = element_blank(),
      axis.line = element_line(size = 0.5),
      axis.ticks =  element_line(size = 0.3),
      axis.text.y = element_text(
        size = 14,
        angle = 0,
        face = "plain"
      ),
      axis.title.x = element_text(
        size = 16
      ),
      #panel.border = element_borderrect(colour = "black", fill=NA, size=5),
      #panel.border = element_border(c("left","bottom")),
      legend.title = element_blank(),
      # switch off the legend title
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.0, "lines"),
      axis.title.y = element_text(size = 12)
    )
  print(XmR_Plot2)
  ggsave(filename = paste0(directory,"/QC_shewhart_plot/","XmR_Plot_polar_",x,".png"), XmR_Plot2,
        width = 42, height = 22, dpi = 300, units = "cm", device='png')
  
}
#View(sel_ISTD)
