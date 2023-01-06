                  

                  #inputData_tto=unique(gsub("(.*)\\D(\\d+$)","\\1",rownames(inputData)))
                  inputData_tto=unique(gsub("(.*)#(.*)\\D(\\d+$)","\\1",rownames(inputData)))
                  inputData_tto=as.data.frame(cbind(Name=unique(inputData_tto),toRemove=NA))
                  
                  fix(inputData_tto)
                  if(length(inputData_tto[!is.na(inputData_tto$toRemove),"Name"])>0){
                    inputData<-inputData[grepl(paste0("^",as.character(inputData_tto[is.na(inputData_tto$toRemove),"Name"]),collapse="|"),rownames(inputData)),]
                  }
                  #go to a temp dir to avoid permission issues
                  inputDataDownload=c()
                  #inputData = t_rsALL.concTOistd.i
                  #eindex = rsALL.samples[,c("index","NAME")]
                  
                  inputData$GRUPO=gsub("(^.*)#(.*)\\D(\\d+$)","\\1",rownames(inputData))
                  inputData$ID=gsub("(^.*)\\D(\\d+$)","\\2",rownames(inputData))
                   #inputData$ID=gsub("(^.*)#(.*)","\\2",rownames(inputData))
                  inputData=cbind(type="SAMPLE",Batch.number=1:nrow(inputData),inputData)
                  
                  colnames(inputData)[3:(ncol(inputData) - 2)] = str_trim(eindex[match(colnames(inputData)[3:(ncol(inputData) -
                                                                                                               2)], eindex$index), "NAME"])                  
                  
                  ### carpeta de descargas
                  directory=get("directory",envir=.GlobalEnv)
                  print("Estas son las carpetas que hay ahora mismo creadas para analisis:")
                  print(matrix(gsub("(.*)\\/(.*)$","\\2",list.dirs(directory))))
                  carpeta_para_figuras=c("ResultadosExcel")
                  names(carpeta_para_figuras) <-c("Introduzca un nombre para la carpeta donde se van a guardar las figuras generadas: ")
                  fix(carpeta_para_figuras)
                  #creamos directorio nuevo
                  dir.create(sprintf(paste0(directory,"/",carpeta_para_figuras)),showWarnings=F);
                  sub_folder<-paste0(directory,"/",carpeta_para_figuras)
                  print(paste("Directorio",sub_folder,"creado"))
                  
                  
                  ## haciendo estadisticas de especie
                  inputDataISTD = summaryRESSTAT(
                    inputData,
                    "GRUPO",
                    colnames(inputData)[-c(1, 2, (ncol(inputData) - 1), ncol(inputData))],
                    varTOcast = "GRUPO",
                    digitos = 3,
                    estat = "media_se",
                    units=paste0("media ",custom_units_to_use)
                  )
                  wide_stats=colnames(inputDataISTD$SummaryStats)[-c(1:3)]
                  inputDataISTD$SummaryStats_wide=as.data.frame(cbind("Sample.Name" = rownames(inputDataISTD$wideData),
                                                                   inputDataISTD$wideData))[,-c(1:2,4)]
                  for(y in wide_stats){
                    aa=reshape2::dcast(inputDataISTD$SummaryStats,GRUPO~Lipid.Name,value.var = y)
                    colnames(aa)[1]=y
                    aa=rbind("-",colnames(aa),aa)
                    colnames(aa)=colnames(inputDataISTD$SummaryStats_wide)
                    inputDataISTD$SummaryStats_wide=rbind(inputDataISTD$SummaryStats_wide,aa)
                  }
                  
                  #View(inputDataISTD$SummaryStats_wide)
                  inputDataISTD$wideData = cbind("Sample.Name" = rownames(inputDataISTD$wideData),
                                                 inputDataISTD$wideData)
                  
                  
                  QC_method=c("ISTD/Cal")
                  names(QC_method)=c("Que analisis estas haciend? ISTD o Cal")
                  fix(QC_method)
                  
                  names(inputDataISTD) = paste0(names(inputDataISTD), "_",QC_method)
                  inputDataDownload = c(inputDataDownload, inputDataISTD)

                 ### haciendo estadisticas de sumas
                   inputDataISTD_class = summaryRESSTAT(
                     inputData,
                     "ID",
                     colnames(inputData)[-c(1, 2, (ncol(inputData) - 1), ncol(inputData))],
                     varTOcast = "GRUPO",
                     digitos = 3,
                     sum = T,
                     estat = "media_se"
                   )
                   colnames(inputDataISTD_class$SummaryStats)[-c(1:3)]=wide_stats
                   inputDataISTD_class$SummaryStats_wide=inputDataISTD_class$wideData
                   for(y in wide_stats){
                     aa=reshape2::dcast(inputDataISTD_class$SummaryStats,GRUPO~Lipid.Name,value.var = y)
                     colnames(aa)[1]=y
                     aa=rbind("-",colnames(aa),aa)
                     colnames(aa)=colnames(inputDataISTD_class$SummaryStats_wide)
                     inputDataISTD_class$SummaryStats_wide=rbind(inputDataISTD_class$SummaryStats_wide,aa)
                    }
                   
                   
                   inputDataISTD_class$wideData = cbind(
                     "Sample.Name" = rownames(inputDataISTD_class$wideData),
                     inputDataISTD_class$wideData
                   )

                   names(inputDataISTD_class) = paste0(names(inputDataISTD_class), "_",QC_method,"_class")
                   inputDataDownload = c(inputDataDownload,
                                         inputDataISTD_class)
                   
     


              ###
                   fs=c()
                   
               for (i in 1:length(inputDataDownload)) {
                 path <- paste(sub_folder,"/",path_sanitize(names(inputDataDownload)[i]), ".csv", sep = "")
                 
                 fs <- c(fs, path)
                 
                 ##### cant be list, converting to character dataframe
                 l <- inputDataDownload[[i]]
                 df = t(data.frame(matrix(
                   unlist(l), ncol = max(lengths(l)), byrow = TRUE
                 )))
                 colnames(df) = colnames(l)
                 rownames(df) = rownames(l)
                 #####
                 
                 write.table(
                   df,
                   file = path,
                   sep = ",",
                   dec = ".",
                   quote = T,
                   fileEncoding  = "UTF-8",
                   row.names = F,
                   col.names = T
                 )
                 
               }
             



