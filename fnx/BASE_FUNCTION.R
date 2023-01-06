##-------------------------------------------------------
##                         FUNCIONES
##--------------------------------------------------------
path_sanitize <- function(filename, replacement = "_") {
  illegal <- "[/\\?<>\\:*|\":]"
  control <- "[[:cntrl:]]"
  reserved <- "^[.]+$"
  windows_reserved <- "^(con|prn|aux|nul|com[0-9]|lpt[0-9])([.].*)?$"
  windows_trailing <- "[. ]+$"
  
  filename <- gsub(illegal, replacement, filename)
  filename <- gsub(control, replacement, filename)
  filename <- gsub(reserved, replacement, filename)
  filename <- gsub(windows_reserved, replacement, filename, ignore.case = TRUE)
  filename <- gsub(windows_trailing, replacement, filename)
  
  # TODO: this substr should really be unicode aware, so it doesn't chop a
  # multibyte code point in half.
  filename <- substr(filename, 1, 255)
  if (replacement == "") {
    return(filename)
  }
  path_sanitize(filename, "")
}

selectOption<-function(x,name="Input",toRemove=F){
  a=as.data.frame(cbind(Input=x,Select=NA),stringAsFactor=F)
  colnames(a)[2]=name
  fix(a)
  a=get("a",envir=.GlobalEnv)
  if(toRemove){
    a=a[is.na(a[,2]),]
  }else {
    a=a[!is.na(a[,2]),]
  }
  return(as.character(a[,1]))
}

## Funciones de escritura de graficos

myPng <- function(x, width=8, height=8, res=300, point=12,dir=getwd()) {
  png(sprintf("%s/%s.png", dir, x), width=width*res, height=height*res, res=res,pointsize=point);
}
myPDF <- function(x, width=7, height=7, pointsize=12, paper="a4",dir=getwd()) {
  pdf(sprintf("%s/%s.pdf", dir, x), width=width, height=height, 
      pointsize=pointsize, paper=paper,useDingbats=TRUE);
}

my.read.csv<-function(x){
  read.csv(x,
           stringsAsFactors = F,
           strip.white = TRUE,
           sep=";",
           na.strings = c("NA","<NA>","#N/A") )
}
my.read.csv2<-function(x){
  read.csv(x,
           stringsAsFactors = F,
           strip.white = TRUE,
           sep=",",
           na.strings = c("NA","<NA>","#N/A") )
}
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
myfun<-function(x){sum(x)}
myfun2<-function(x) sum(x[(x>0)&!is.na(x)]);#si se quiere evitar los NA
f.Norm2<-function(x){x/sum(x[(x>0)&!is.na(x)])};
f.Norm<-function(x){x/sum(x)};
f.Norm_sc<-function(x){scale(x,scale=T)}; #autoscaling
f.Norm_mc<-function(x){scale(x,scale=F)}; #mean-centered scaling
my.mean<-function(x){mean(x[(x>0)&!is.na(x)])};
my.sd<-function(x){sd(x[(x>0)&!is.na(x)])};
my.max<-function(x){max(x[(x>0)&!is.na(x)])};
my.min<-function(x){min(x[(x>0)&!is.na(x)])};
my.median<-function(x){median(x[(x>0)&!is.na(x)])};
my.median2<-function(x){median(x[!is.na(x)])};
my.mean.rep<-function(x){
  rep(signif(mean(x[(x>0)&!is.na(x)]),digits=6)
             ,times=length(x))
}
my.N<-function(x){length(x[(x>0)&!is.na(x)])}; #count the number of occurrences

#Imputation of zero values to the value passed in fuction it can be and array given
my.imputation<-function (x,value,na.rm=FALSE){
  if (na.rm){
    for (i in c(1:nrow(x))){    
      for (j in c(1:ncol(x))){
        
        ifelse(!is.na(x[i,j]) & (x[i,j]==0),x[i,j]<-value[j],1)
        #to ONLY change 0 values let missing values NA unchanged, na.rm=FALSE default
        
      }
    }}
  else{
    for (i in c(1:nrow(x))){    
      for (j in c(1:ncol(x))){
        
        ifelse(is.na(x[i,j]) | (x[i,j]==0),x[i,j]<-value[j],1)
        #to change both 0 and missing values NA by criteria, na.rm=TRUE
        
      }
    } 
  }
  return(x)    
}

##graficos de visualizacion de normalizaciones
my.graphNorm<-function(h){
  
  bp<-boxplot(h,axes=F,frame.plot=T,cex.axis=1.0);
  stripchart(h,vertical=TRUE,pch=1,add=T,cex=0.2);
  axis(1, at = c(1:ncol(h)), labels =bp$names, 
       lwd=.5, #longitud de los tick 
       lwd.ticks=1.0,
       las=2,  #orientacion
       cex.axis=.5, #tama??o de los label
       tck=0.02,
       mgp=c(0, 0.02, 0)
  );
} 

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <-
  function(data = NULL,
           measurevar,
           groupvars = NULL,
           na.rm = FALSE,
           conf.interval = .95,
           .drop = TRUE,
           notification = F) {
    require(plyr)
    
    #save objects in current environment
    # save(list = ls(), file = "shiny_envse.Rdata", envir = environment())
    # load("shiny_envse.Rdata")
    data = data[, c(groupvars, measurevar)]
    colnames(data) = c(groupvars, "klx")
    datac = data %>%
      group_by_at(vars(groupvars), .drop = T) %>% summarise(
        N = sum(!is.na(klx)),
        mean = my.mean(klx),
        sd = my.sd(klx, notification = notification),
        median = my.median(klx),
        suma = myfun2(klx),
        cv = 100 * sd / mean,
        se = sd / sqrt(N),
        ci = sd / sqrt(N) * qt(conf.interval / 2 + .5, N - 1)
      )
    colnames(datac)[which(colnames(datac) %in% "mean")] <- measurevar
    return(as.data.frame(datac))
  }

my.max <- function(x)
  ifelse(!all(is.na(x)), max(x, na.rm = T), NA)
myfun <- function(x) {
  sum(x, na.rm = T)
}
myfun2 <-
  function(x)
    sum(x[(x > 0) & !is.na(x)])
#si se quiere evitar los NA
f.Norm2 <- function(x) {
  x / sum(x[(x > 0) & !is.na(x)])
}

f.Norm <- function(x) {
  x / sum(x)
}

f.Norm_sc <- function(x) {
  scale(x, scale = T)
}
#autoscaling
f.Norm_mc <- function(x) {
  scale(x, scale = F)
}
#mean-centered scaling
my.mean <- function(x) {
  mean(x[(x > 0) & !is.na(x)])
}

my.sd <- function(x, notification = F) {
  if (sum(is.na(x)) == length(x)) {
    return(0)
  }
  sd = sd(x[(x > 0) & !is.na(x)])
  if (is.na(sd) | sd == 0) {
    sd = x * 0.15
    if (notification == T) {
      showNotification(
        "Too few replicates. RSD set to 15% by default.",
        type = "warning",
        duration = 5
      )
      
    }
  }
  return(sd)
}

my.max <- function(x) {
  max(x[(x > 0) & !is.na(x)])
}

my.min <- function(x) {
  min(x[(x > 0) & !is.na(x)])
}

my.median <- function(x) {
  median(x[(x > 0) & !is.na(x)])
}

my.median2 <- function(x) {
  median(x[!is.na(x)])
}

my.mean.rep <- function(x) {
  rep(signif(mean(x[(x > 0) & !is.na(x)]), digits = 2)
      , times = length(x))
}

# Correlation matrix with p-values. See http://goo.gl/nahmV for documentation of this function
cor.prob <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}

## Use this to dump the cor.prob output to a 4 column matrix
## with row/column indices, correlation, and p-value.
## See StackOverflow question: http://goo.gl/fCUcQ
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}

MyInclusionCriteria <- function (xx,umbral=4,nas=T,factor=1) {
  # nas=T let NAs without imputation impute only zero values
  ##recibe un vector y lo transforma seg??n el umbral (n?? observaciones minimas para hacer imputacion)
  #t_matrix<-t_rsALL.plasma[,-c(1:2)];
  #factor, aplica una correccion al nivel minimo para salvaguardar diferencias extremas en algunos casos
  zero_val<-data.frame(apply(xx,2,my.N)); #vector con datos que NO son 0 o NA
  index<-xx[,which(zero_val>=umbral)]
  min_val<-as.vector(apply(xx[,which(zero_val>=umbral)],2,my.min)/factor); #vector con datos del valor minimo x columna
  xx[,which(zero_val>=umbral)]<-my.imputation(index,min_val,na.rm=nas);  #Transforming zero values to min values na.rm=TRUE to let NAs
  ### el resto son NAs, pero no pierdo el valor de los que hay. Esto no se si es bueno
  xx[xx==0]<-NA; 
  return(xx);
}


###Function to trasnform by ISTD concentration####

MyMatrix_to_istd <- function (matrix1,matrix2,multipliers=NULL) {
  ## get matrix of areas and transform to concentrations
  ## matrix1 is the global matrix
  ## matrix2 is the matrix of calibrators/ISTDs
  ## multipliers is the conc factors, if not given then is 1 by defect
  matrix1_norm<-NULL
  matrix_aux<-NULL
  factors<-levels(matrix2[,1])
  for (j in factors){
    aux<-matrix2[which(matrix2[,1]==j),-c(1)]
    if (nrow(aux)==0){aux<-rep(1,length(aux))} #handle if not ISTD
    matrix1_norm<-matrix1[which(matrix1[,1]==j),-c(1)]
    matrix1_norm<-sweep(as.matrix(matrix1_norm),2,as.matrix(aux),FUN="/",check.margin=T)
    multiplier<-c(1)
    if(!is.null(multipliers)){
      multiplier<-multipliers[which(multipliers[,1]==j),-c(1)]}
    #if (nrow(multiplier)==0){aux<-rep(1,length(aux))} #handle if not ISTD    
    matrix1_norm<-matrix1_norm*as.numeric(multiplier)
    rownames(matrix1_norm)<-rep(j,length(matrix1_norm[,1]))
    #matrix_aux<-matrix1_norm
    #print(matrix1_norm)
    matrix_aux<-rbind(matrix_aux,matrix1_norm)
  } 
  return(matrix_aux)
}

hboutlier <- function(x,r=3){
  x[!is.finite(x)] <- 0.00000000000000000000001 #para mantener el formato de output de Detectoutliers, sustituyo Na
  x <- x[is.finite(x)]
  stopifnot(
    length(x) > 0
    , all(x>0) )
  xref <- median(x)
  if (xref <= sqrt(.Machine$double.eps))
    warning("Reference value close to zero: results may be inaccurate")
  pmax(x/xref, xref/x) > r
}

###Function to detect ouliers####
DetectOutliers<-function(x_out,r=2){
  ###x_out is a data.frame of tipe t_rsALL
  z<-list();
  for (i in c(1:ncol(x_out))){
    #z[[a]]<-hboutlier(x_out[,i],r=2) # r=4 times upper whisker with function outliers
    #outliers<-round(x_out[,i],digits=4)==round(boxplot.stats(x_out[,i],coef=r)$out,digits=4)
    #outliers<-x_out[,i]==boxplot.stats(x_out[,i],coef=r)$out    
    outliers<- x_out[,i] %in% boxplot.stats(x_out[,i],coef=r)$out 
    if (length(outliers)>0){
      z[[i]]<-outliers 
    } else { z[[i]]<-as.logical(rep("FALSE",times=length(x_out[,i]))) }
    #boxplot(x_out[,i])
    #x_out[which(hboutlier(x_out[,i],r=2)==T),i]
  } 
  z<-matrix(unlist(z), ncol = length(z), byrow = F)
  return(z)
}


# Function for plotting colors side-by-side
pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}



To_Numeric_Matrix<-function(s){
  coln<-colnames(s)
  rown<-rownames(s)
  ## get char matrix with especial characters and convert to numeric sanitizing extrange chars
  xx<-list();
  for (i in c(1:length(s))){
    xx[[i]]<-as.numeric(sub("<","",as.character(s[,i])))
  }
  xx<-matrix(unlist(xx), ncol = length(s), byrow = F)
  colnames(xx)<-coln
  rownames(xx)<-rown
  return(xx)
  
}

isotopic_correction_type1 <-
  function(rsALL, db) {
    # save(rsALL,db,file="iso_test.rda")
    # load("iso_test.rda")
    #db = rsISTD[, c("NAME", "NCarbon")]
    rsALLbkp = rsALL
    rsALL$NAME = gsub("\\[.*", "", rsALL$NAME)
    rsALL$NAME = sub("d\\d{1,2}", "", rsALL$NAME)
    
    db$NAME = str_trim(db$NAME)
    db=db[which(str_trim(db$NAME)%in% str_trim(rsALL$NAME) ),]
    classes = which(str_trim(rsALL$NAME) %in% str_trim(db$NAME))
    matcheNames = match(str_trim(db$NAME) ,str_trim(rsALL$NAME))
    db$classes = rsALL[matcheNames, "Region"]
    
    
    # rsALL$NAME = gsub("(.*\\(.*)(\\(.*\\))(.*\\))", "\\1\\3", rsALL$NAME)
    # db$NAME = gsub("(.*\\(.*)(\\(.*\\))(.*\\))", "\\1\\3", db$NAME)
    
    withParentesis = grepl("\\d{1,2}\\:.*\\[.*", rsALL$NAME)
    if (sum(withParentesis) > 0) {
      rsALL[withParentesis, "NAME"] = gsub("\\[.*", "", rsALL[withParentesis, "NAME"])
      rsALL[!withParentesis, "NAME"] = gsub("\\[|\\]", " ", rsALL[!withParentesis, "NAME"])
    } else {
      rsALL$NAME = gsub("\\[|\\]", "", rsALL$NAME)
      
    }
    db[!is.na(db$classes), "NAME"] = rsALL[matcheNames[!is.na(matcheNames)], "NAME"]
    
    
    
    
    carbon = gsub("\\.\\d|:\\d|[[:alpha:]]+|\\s$|\\s^|(\\d+)(?!\\:|\\d)", "", rsALL$NAME,perl=T)
    carbon = gsub("\\D|\\s", " ", carbon)
    carbon = gsub("(?<=\\d)(\\s)(?=\\d)", "+", carbon,perl=T)
    carbon = gsub("\\s", "", carbon)
    carbon = gsub("\\+", " ", carbon)
    carbon = strsplit(carbon, " ")
    carbon = sapply(carbon, function(x)
      as.numeric(as.character(x)))
    carbon = sapply(carbon, sum)
    
    isotope_typeI <-
      as.data.frame(cbind(
        NAME = as.character(rsALL$NAME),
        class = as.character(paste0(rsALL$Region, "=")),
        carbon = as.numeric(carbon)
      ),stringsAsFactors =F)
    isotope_typeI[, c(2)] <- sub("cyl", "", isotope_typeI[, c(2)])
    isotope_typeI <- data.frame(isotope_typeI)
    isotope_typeI[which(str_trim(isotope_typeI$NAME) %in% str_trim(db$NAME)), "class"] =
      gsub("\\=", "-IS", isotope_typeI[which(str_trim(isotope_typeI$NAME) %in% str_trim(db$NAME)), "class"])
    classes = classes[!is.na(classes)]
    
    classes = as.data.frame(cbind(
      class = as.character(rsALL[classes, c("Region")]),
      NAME = as.character(rsALL[classes, c("NAME")]),
      NCcalculated = as.numeric(carbon[classes])
    ),stringsAsFactors =F)
    classes$NAME = str_trim(classes$NAME)
    db = db[!is.na(db$classes), ]
    db$NAME = str_trim(db$NAME)
    db = merge(classes, db, id = "NAME",)
    db$NCarbon = as.numeric(db$NCarbon) - as.numeric(db$NCcalculated)
    db$NCarbon = ifelse(db$NCarbon < 0, 0, db$NCarbon)
    
    for (i in unique(db[, "class"])) {
      #i="ACer"
      Ncarbon <-
        as.numeric(as.character(isotope_typeI[grepl(paste0(i, "\\=|", i, "-IS"), isotope_typeI[, "class"]), "carbon"])) +
        as.numeric(db[grepl(i, db$class), "NCarbon"])
      #which(isotope_typeI[grepl(paste0(i, "\\=|", i, "-IS"), isotope_typeI[, "class"]), "NAME"] %in% db$NAME)
      z_factor <-
        sapply(Ncarbon, function(x) {
          1 + 0.0109 * x + 0.0109 ^ 2 * x * (x - 1) / 2
        }) ## calular factor Z type I correccion isotopica
      isotope_typeI[grepl(paste0(i, "\\=|", i, "-IS"), isotope_typeI[, "class"], ignore.case = T), "Ncarbon"] <-
        Ncarbon
      isotope_typeI[grepl(paste0(i, "\\=|", i, "-IS"), isotope_typeI[, "class"], ignore.case = T), "Z"] <-
        z_factor  #factor corrector areas
      ##correccion isotopica para concentraciones
      isotope_typeI[grepl(paste0(i, "\\=|", i, "-IS"), isotope_typeI[, "class"]), "Z1"] <-
        isotope_typeI[grepl(paste0(i, "\\=|", i, "-IS"), isotope_typeI[, "class"]), "Z"] / isotope_typeI[grepl(paste0(i, "-IS"), isotope_typeI[, "class"], ignore.case = T), "Z"]
      
    }
    isotope_typeI[is.na(isotope_typeI[, "Z1"]), "Z1"] = 1
    
    # isotope_typeI$lipid_abbrev=paste0(isotope_typeI$class,"(",carbon,":",insaturation)
    # library(jsonlite)
    # json=sapply(isotope_typeI$lipid_abbrev, function(x) fromJSON(paste0("https://www.lipidmaps.org/rest/compound/abbrev/",x,"/formula/json"))$Row1$formula)
    #
    #
    
    
    rsALL.isotopic <- NULL
    if (all(isotope_typeI$NAME == rsALL$NAME)) {
      rsALL.isotopic <-
        sweep(rsALL[, -c(1:which(colnames(rsALL) %in% c("ERROR")))], 1, isotope_typeI$Z, FUN =
                "*")
      rsALL.isotopic <-
        cbind(rsALL[, c(1:which(colnames(rsALL) %in% c("ERROR")))],rsALL.isotopic)
    } else {
      stop("NO SE REALIZ? LA CORRECCION ISOTOPICA")
    }
    
    rsALL.isotopic = rsALL.isotopic[match(as.numeric(rownames(rsALL.isotopic)), as.numeric(rownames(rsALLbkp))), ]
    rsALL.isotopic$NAME = rsALLbkp$NAME
    
    return(rsALL.isotopic)
  }
options(scipen=999) #set scientific notation off
#options(scipen=0) #set scientific notation on
options(digits=5) #set the number of digits default is 7