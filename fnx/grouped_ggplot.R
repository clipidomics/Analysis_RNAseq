library(ggpubr)
library(rstatix)
# Transform `dose` into factor variable
# dfe <- ToothGrowth
# dfe$dose <- as.factor(df$dose)
# head(dfe, 3)

###transform to long format
dfc.long<-species_barplot_data
head(dfc.long, 3)
custom_order<-c("NORMAL","NAFL","NASH","NASHFSIG")
escala=c("white",c(brewer.pal(9,"Blues")))[c(1,3,5,7,9,10)]

dfc.long$tto<-factor(dfc.long$tto,levels=custom_order,ordered=T)
referencia<-"NORMAL" ##set to NA si no se quiere 
step.upper<-1.15 # salto de eje-y para que queden bien las barras de comparaciones
#!is.na(referencia)
species_plot_list<-list()
stat_test_list<-list()
for (cls in region){
  
      #cls<-"dhHexCer"
      #creamos un subset con regiones indicadas arriba
      sub_fusion<-subset(dfc.long, dfc.long$TAG %in% cls & dfc.long$tto %in% custom_order)
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
      if(!is.na(referencia)){
        stat.test <- stat.test %>%
        add_xy_position(fun = "mean_sd", x = "m_NAME", dodge = 0.8, ref.group=referencia,step.increase = 0.05) 
      }else{
       stat.test <- stat.test %>%
        add_xy_position(fun = "mean_sd", x = "m_NAME", dodge = 0.8,step.increase = 0.05) 
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
      # +geom_dotplot(
      #   binaxis = 'y',
      #   stackdir = 'center',
      #   stackratio = 0.8,
      #   dotsize = .4,
      #   alpha = 0.8,
      #   aes(x=m_NAME,y = value.y, fill = tto),
      #   #aes(x=m_NAME,y = value.y),
      #   show.legend = FALSE,
      #   position = position_dodge(.8)
      # )
      # 
      bp<-bp + #facet_wrap( ~ variable, scales="free",ncol=3)+ ##Coment to separate classes
        xlab("") +
        #ylab(paste("% of total", sub_fusion[1,3], sep=" "))+ #ojo con unidades
        ylab("") +
        #geom_hline(aes(yintercept = 0),
        #          colour = "black",
        #         linetype = "solid") +
        theme_classic(base_line_size = 1) +
        labs(title = cls) +
        scale_y_continuous(expand=c(0,0),limits=c(0,max(sub_fusion$value.y)*step.upper))+
        #scale_y_continuous(expand=c(0,0))+
        #scale_x_continuous(expand = expansion(mult = c(0, .2))) +
        #guides( fill = FALSE)+#labs(subtitle = x)+
        #guides(shape=F)+
        theme(
          #axis.text.x = element_blank(),
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
      
      species_plot_list[[cls]]<-bp
      stat_test_list[[cls]]<-stat.test
      print("?Barplot de especies hecho!")
      
      #print(bp)
      
      #pushViewport(viewport(y=0.95,x=0.85,height=.05))
      #title=unlist(strsplit(directory,"/"))[length(unlist(strsplit(directory,"/")))-1]
      #grid.table(title)
      
      #dev.off()
      }
      
}

#species_plot_list

ggsave(filename = paste0(directory,"/fig_molec_PLASMA_annotated_stats",".png"),
       annotate_figure(
         ggarrange(species_plot_list$dhCer, species_plot_list$Cer,
                   species_plot_list$dhSM,species_plot_list$SM,
                   species_plot_list$dhHexCer, species_plot_list$HexCer,
                   #labels = c("A", "B", "C","D","E","F"),
                   ncol = 2, nrow = 3) #,
         #left = text_grob("nmol/mL", color = "black", face = "plain", size = 16,rot=90, vjust=3.5)
       ),
       width = 22, height = 22, dpi = 600, units = "cm", device='png')

ggsave(filename = paste0(directory,"/fig_molec_PLASMA_neutral_stats",".png"),
       annotate_figure(
         ggarrange(species_plot_list$CE, species_plot_list$FC,
                   species_plot_list$TG,species_plot_list$ACer,
                   species_plot_list$PC, species_plot_list$PE,
                   #labels = c("A", "B", "C","D","E","F"),
                   ncol = 2, nrow = 3) #,
         #left = text_grob("nmol/mL", color = "black", face = "plain", size = 16,rot=90, vjust=3.5)
       ),
       width = 22, height = 22, dpi = 600, units = "cm", device='png')
print("?Barplot de especies hecho!")



