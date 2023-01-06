#inputData=t_rsALL.concTOistd.i

######
## Heatmap test
######
inputData_tto=unique(gsub("(.*)\\D(\\d+$)","\\1",rownames(inputData)))
inputData_tto=as.data.frame(cbind(Name=unique(inputData_tto),toRemove=NA))

QC_method=c("ISTD/Cal")
names(QC_method)=c("Que analisis estas haciend? ISTD o Cal")
fix(QC_method)

heatmapType=c("Suma/Especie")
names(heatmapType)=c("Que tipo de grafico quieres ver? Suma de clases o de especies individuales")
fix(heatmapType)


fix(inputData_tto)
if(length(inputData_tto[!is.na(inputData_tto$toRemove),"Name"])>0){
  inputData<-inputData[!grepl(paste0("^",as.character(inputData_tto[!is.na(inputData_tto$toRemove),"Name"]),".*",collapse="|"),rownames(inputData)),]
}
#### quitar classes
inputData_class=unique(gsub("(.*)\\D(\\d+$)","\\1",colnames(inputData)))
inputData_class=as.data.frame(cbind(class=unique(inputData_class),toRemove=NA))

fix(inputData_class)
if(length(inputData_class[!is.na(inputData_class$toRemove),"toRemove"])>0){
  
  inputData<-inputData[!grepl(paste0("^",as.character(inputData_class[!is.na(inputData_class$toRemove),"class"]),".*",collapse = "|"),colnames(inputData))]
}


#### secuencia para generar un dfc.sum
dfc = inputData
dfc$tto = gsub("\\..*", "", rownames(dfc))
dfc = reshape2::melt(dfc, id = c("tto"))

if(heatmapType=="Suma"){
  dfc$variable = gsub("(\\D+)_(\\d+$)", "\\1", dfc$variable)
}

dfc = ddply(dfc, .(variable, tto), colwise(myfun2, c("value")))
dfc$ID = gsub(".*(\\D+)(\\d+$)", "\\2", dfc$tto)
dfc$tto = gsub("(\\D+)(\\d+$)", "\\1", dfc$tto)
dfc.sum = summarySE(dfc,
                    measurevar = "value",
                    groupvars = c("tto", "variable"))



head(dfc.sum)
dfc.sum[(dfc.sum$cv < 25), "cv"] = NA
dfc.sum_normalize_by=selectOption(as.character(unique(dfc.sum$tto)),"Sample_GROUP_to_normalize_by")
dfc.sum=dfc.sum%>%group_by(variable)%>%summarise(tto=tto,cv=cv,value=value/my.mean(value[grepl(paste0("^",dfc.sum_normalize_by),tto)]))
dfc.sum[, "value"] = round((dfc.sum[, "value"]-1)/0.2,digits=0)*20




bb = reshape2::dcast(dfc.sum[!grepl("CAL|QC", dfc.sum$tto), c("tto", "variable", "cv")], tto ~ variable, value.var = "cv")
rownames(bb) = gsub("(^.*_)(\\d+$)", "\\2", bb$tto)

aa = reshape2::dcast(dfc.sum[!grepl("CAL|QC", dfc.sum$tto), c("tto", "variable", "value")], tto ~ variable, value.var = "value")
rownames(aa) = gsub("(^.*_)(\\d+$)", "\\2", aa$tto)



normalized_values=paste0(min(dfc.sum[, "value"]),",",max(dfc.sum[, "value"]))
names(normalized_values) = paste0(
  "Escala de valores: ",
  paste0(as.character(sort(
    as.numeric(unlist(unique((
      dfc.sum[, "value"]
    ))))
  )), collapse = ", "),
  ". Indica el rango de interes"
)
fix(normalized_values)
normalized_values=unlist(str_split(normalized_values,","))


gradient_col <- ggplot2::scale_fill_gradient2(
  low = "#000000",
  high = "#008000",
  mid = "#dbdcdd" ,
  na.value = "#FFFFFF", #color por encima de la escala
  midpoint = 0,
  limits = c(as.numeric(normalized_values[1]),as.numeric(normalized_values[2])) ## valores por encima se veran como na.value
)

clases=gsub("_\\d+$","",colnames(aa))[-1]
if(heatmapType=="Especie"){
  colnames(aa)=c("tto",eindex[match(colnames(aa)[-1],eindex$index),"NAME"])
}
##cosas para poder guardar la imagen
library(heatmaply)
#webshot::install_phantomjs()
p=heatmaply((aa),
          Colv = F,
          Rowv = T,
          seriate = "GW",
          #dendrogram = F,
          grid_gap = 1,
          scale_fill_gradient_fun =(gradient_col),
          col_side_colors =clases,
          file = paste0("heatmaply_plot_",QC_method,"_",heatmapType,".html"),
          label_names = c("ID", "LIPID", "% respect control group"),
          xlab = "The numbers correspond to the CV over 25%.",
          cellnote = round(bb[, -1], digits = 0),
          cellnote_textposition = "bottom center",
          cellnote_size = 10,
          colorbar_yanchor = "bottom",
          column_text_angle =90,
          showticklabels=c(F, T)
          #colors = palete
)
print(p)
print("Heatmap saved in root directory as 'heatmaply_plot.html<'")
####
## heatmap ends
####