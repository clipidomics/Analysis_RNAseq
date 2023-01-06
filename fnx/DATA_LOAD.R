### necesita dos variables para indicar que archivos coger
data_files


## empieza aqui
  rm(input)
  ### Fusionamos todos los archivos con resultados
  for(file in data_files){
     #file<-"/Users/opastor/OneDrive - Madrid Digital/YALE/PROYECTO ASOS_DEGS1_2/5_RT_PCR/Analysis RT_PCR_06102022/Degs1 y Degs2 y 18s liverOscar -  Quantification Summary_06102022.csv"
    input_temp<- data.frame(lapply(file,my.read.csv))
    fix(input_temp)
    ###HACER 
     if(sum(colnames(input_temp)=="Well")>0){
       colnames(input_temp)[1] <- "PCR"
       # transformar datos a formato conocido
       input_temp<-cbind(input_temp,GRUPO=gsub("(.*)\\_(.*)","\\2",input_temp$Sample))
       input_temp<-cbind(input_temp,Muestra=gsub("(.*)\\_(.*)","\\1",input_temp$Sample))
       qpcr_order<-c("PCR","GRUPO","Muestra","Sample","Content","Target","Cq")
       input_temp<-input_temp[qpcr_order]
       input_temp<-input_temp[which(!input_temp$GRUPO==""),]
       #input_temp_old<-input_temp
       input_temp<-input_temp%>%group_by(PCR,Sample,Target)%>%summarise(PCR=unique(PCR),
                                                                 GRUPO=unique(GRUPO),
                                                                 Muestra=unique(Muestra),
                                                                 Sample=unique(Sample),
                                                                 Content=unique(Content),
                                                                 Target= unique(Target),
                                                                 Cq=my.mean(Cq)
       )
       
       input_temp<-input_temp[qpcr_order]  
     }
    
      ## seleccionamos datos base 
      input_temp=input_temp[,c(1:7)]
      colnames(input_temp)=c("Analysis","GRUPO","Muestra","Sample.Name","Type","Gen","Mean.Cp")
      input_temp[1:6]=lapply(input_temp[1:6],as.character)
      input_temp[7:ncol(input_temp)]=lapply(input_temp[7:ncol(input_temp)],as.numeric)
      
      
      
      ## creamos la media de los genes HK
      Mean.HK=input_temp[which(input_temp$Type=="Housekeeping"),]%>%group_by(Muestra,Analysis)%>%summarise(Analysis=unique(Analysis),
                                                                                                   GRUPO=unique(GRUPO),
                                                                                                   Muestra=unique(Muestra),
                                                                                                   Sample.Name=unique(Sample.Name),
                                                                                                   Type=unique(Type),
                                                                                                   Gen="Mean.HK",
                                                                                                   Mean.Cp=my.mean(Mean.Cp)
                                                                                                   )
      n = colnames(input_temp)
      #input_temp=join(x, y, by = NULL, type = "left", match = "all")
      input_temp=rbind(data.frame(input_temp[n]),data.frame(Mean.HK[n]))
      rm(Mean.HK)
      
      
      ## buscamos HK que hay
      HK_found=unique(input_temp[which(input_temp$Type=="Housekeeping"),"Gen"])
      
      ## hacemos el calculo 
      for(hk in HK_found){
        input_temp=input_temp%>%group_by(Muestra,Analysis)%>%mutate(value=(1/(2^(Mean.Cp-Mean.Cp[which(Gen==hk)[1]]))))
        colnames(input_temp)[ncol(input_temp)]=hk
      }
      
      
      #guardamos los cambios
      if(!exists("input")){
        input <- input_temp
      } else {
        
        ### corregimos diferencias en los nombres
        input_temp_genes=input$Gen[match(tolower(input_temp$Gen),tolower(input$Gen))]
        ## los que son NAs los conservamos
        input_temp_genes[is.na(input_temp_genes)]=input_temp$Gen[is.na(input_temp_genes)]
        input_temp$Gen=input_temp_genes
        
        input <- rbind(input,input_temp)
      }
      rm(input_temp)
      
      
  }
  
### quitamos muestras vacias
  input=input[!is.na(input$Mean.Cp),]

