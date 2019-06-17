# Mehrnoosh Oghbaie
# 06/12/2019
# Repository for all the functions
# This file includes most of the general functions that are going to be used in this project

######################################################################################
# Either download or install the required library from CRAN or bioconductor
#######################################################################################

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  #if (length(bioconductor.packages) > 0) {
  #  source("http://bioconductor.org/biocLite.R")
  #}
  
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      BiocManager::install(p,version = "3.8") 
    }
    library(p,lib.loc="~/R/win-library/3.5",character.only=T)
  }
  
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p) 	
    }	
    #library(p, lib.loc="~/R/win-library/3.5",character.only=T) 
    library(p, lib.loc="~/R/win-library/3.5",character.only=T)
  }
}



###################################################################################
### Imputation functions
##################################################################################

### function : finding column with minimum zeros
count_zeros <- function(x){
  return(sum(is.na(x)))
}

min_col <- function(x){
  return(names(x[x==min(x)]))
}
### function: generate uniform distribution for given number, mean and std
na_zeros_impute <- function(na_zeros , mu,sd){
  return(matrix(runif(sum(na_zeros), mu-3*sd, mu-2*sd),ncol=1))
}

perseus_zeros_impute <- function(na_zeros , mu,sd){
  return(matrix(rnorm(sum(na_zeros), mean = mu-1.8*sd, sd = 0.3*sd),ncol=1))
}

#####################################################################################################
### Imputation method 6_2

#Method_2:
#  For each protein, which has zeros in each replica, replica with the least number of zeros is chosen.
#  For it mean and standard deviation of non-zero proteins intensities are counted. 
#  New intensity for this protein in this replica is sampled from uniform distribution with parameters: 
#  start = mu - 3*sd, end = mu - 2*sd. Then other replicas for this protein are imputed as in 6th step
#Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

#Method_6:
# Impute values for proteins, which have zero intensities only in some replicates:
#  Build distribution of deltas for all non zero proteins, where 
# delta =(Intrep1-Intrep2) mean(Intrep1,Intrep2)
# Calculate mudelta , sddelta
# Calculate new delta and new Intensity:
# deltanew=rnorm(mu=mudelta, sd=sddelta*sqrt(2)*mean(correlations))
# Inew=mean(Intother)*abs(1+deltanew)



### Method2
impute_all_zeros_2 <- function(data,con){
  ## condition: two condition in each table
    colnames <- colnames(data)[grepl(con,colnames(data))]
    na_zeros <- rowSums(data[,colnames],na.rm=T)==0
    min_zero_col <- min_col(apply(data[,colnames],2,count_zeros))
    sample <- unname(unlist(data[,min_zero_col]))
    mu <- mean(sample, na.rm=T)
    sd <- sd(sample, na.rm=T)
    data[na_zeros,colnames] <- matrix(na_zeros_impute(na_zeros* length(colnames) , mu,sd),ncol= length(colnames) , nrow = sum(na_zeros))
    return(data)
    }
  

### Funtion that Imputate partial rows
### Method6

impute_partial_zero_6_2 <- function(data,condition){
  for(con in condition){
    data <- impute_all_zeros_2(data,con)
    colnames <- colnames(data)[grepl(con,colnames(data))]
    count_na <- apply(data[,colnames],1, count_zeros)
    y <- data[count_na>0,colnames]
    print(paste("There are ",dim(y)[1]," records missing. \n And ", (dim(data)[1]-dim(y)[1]), " records with full values in ",condition))
    if(length(colnames)>2){
      for(i in 1:dim(y)[1]){
        col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
        for(j in col_NAs){
          col_NA <- j
          col_select <- names(y[i,colnames(y)!=col_NA])
          sample <- y[complete.cases(y[,col_select]),col_select]
          delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
          mu <- mean(unlist(unname(delta)), na.rm=T)
          std <- sd(unlist(unname(delta)),na.rm=T)
          cor <- cor(data[count_na==0,colnames])
          mean_cor <- mean(cor[col_NA,col_select])
          deltanew <- rnorm(1,mu, std*sqrt(2)/mean_cor)
          y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
        }
      }
    }else{
      na_zeros <- length(y[is.na(y)])
      mu <- mean(y[!is.na(y)])
      sd <- sd(y[!is.na(y)])
      y[is.na(y)] <- na_zeros_impute(na_zeros , mu,sd)
    }
    
    data[count_na>0,colnames] <- y
    
  }
  return(data)
} 


#Method_1:
# For each replica mean and standard deviation of non-zero proteins intensities are counted. 
# New intensity for each missing value in each replica is sampled from uniform distribution with parameters:
#  start = mu - 3*sd, end = mu - 2*sd
# Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

#Method_7:
# For outliers (proteins, which have zero values in all but one replica) this method implies one of the methods for imputation of all zero replicas.
# For other proteins uses method6

### function: mutate all NA values for specific condition( control or case)
### Method1
impute_all_zeros_1 <- function(data,condition){
  if(length(condition)>1){
    condition <- paste0(condition, collapse="|")
  }
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  less1zero <- rowSums(data[,colnames],na.rm=T)<=1
  for(i in 1:length(colnames)){
    zero_col <- is.na(data[,colnames[i]])
    sample <- data[!less1zero,colnames[i]]
    sample <- sample[complete.cases(sample)]
    print(paste("There are ", length(sample), " complete records in ", colnames[i], ".\n There are ",sum(unname(less1zero&zero_col)), "zero records to fill"))
    #OutVals = boxplot(sample)$out
    #sample <- sample[!sample %in% OutVals]
    print(paste( length(sample), "records are used for imputing zero values in ", colnames[i]))
    mu <- mean(sample, na.rm=T)
    sd <- sd(sample, na.rm=T)
    data[less1zero&zero_col,colnames[i]] <- na_zeros_impute(less1zero&zero_col, mu,sd)
    #data[less1zero&zero_col,colnames[i]] <- perseus_zeros_impute(less1zero&zero_col, mu,sd)
  }
  return(data)
}


### Method7_1
impute_partial_zero_7_1 <- function(data=list_experiment$Colon,condition="Tumor"){
  print(names(data))
  if(length(condition)>1){
    condition <- paste0(condition, collapse="|")
  }
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  count_na_initial <- apply(data[,colnames],1, count_zeros)
  data <- impute_all_zeros_1(data,condition)
  count_na <- apply(data[,colnames],1, count_zeros)
  y <- data[count_na>0,colnames]
  print(paste("There are ",dim(y)[1]," records missing. \n And ", (dim(data)[1]-dim(y)[1]), " records with full values in ",condition))
  for(i in 1:dim(y)[1]){
    col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
    for(j in col_NAs){
      col_NA <- j
      col_select <- names(y[i,colnames(y)!=col_NA])
      sample <- y[complete.cases(y[,col_select]),col_select]
      delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
      mu <- mean(unlist(unname(delta)), na.rm=T)
      std <- sd(unlist(unname(delta)),na.rm=T)
      cor <- cor(data[count_na==0,colnames])
      mean_cor <- mean(cor[col_NA,col_select])
      deltanew <- rnorm(1,mu, std*sqrt(2)/mean_cor)
      y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
    }
  }
  data[count_na>0,colnames] <- y
  return(data)
} 

### Perseus imputation


impute_Perseus <- function(data,condition){
  if(length(condition)>1){
    condition <- paste0(condition, collapse="|")
  }
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  for(i in 1:length(colnames)){
    zero_col <- is.na(data[,colnames[i]])
    sample <- data[,colnames[i]]
    sample <- sample[complete.cases(sample)]
    print(paste("There are ", length(sample), " complete records in ", colnames[i], ".\n There are ",sum(unname(zero_col)), "zero records to fill"))
    #OutVals = boxplot(sample)$out
    #sample <- sample[!sample %in% OutVals]
    print(paste(length(sample), "records are used for imputing zero values in ", colnames[i]))
    mu <- mean(sample, na.rm=T)
    sd <- sd(sample, na.rm=T)
    data[zero_col,colnames[i]] <- perseus_zeros_impute(zero_col, mu,sd)
  }
  return(data)
}


###################################################################################
### Anova functions
##################################################################################
#### Volcano plot

draw_volcanoplot <- function(data, condition){
  outliers <- PA[["experimentPreImputateSignificant"]][[condition]]$Case_non_zero>1
  ds <- data[[condition]]
  ds[["outliers"]] <- outliers
  ds <- ds[complete.cases(ds),]
  ds$uniprotID <- unlist(lapply(ds$uniprotID, function(x) strsplit(x,";")[[1]][1]))
  ds$uniprotID <- unlist(lapply(ds$uniprotID, function(x) strsplit(x,"-")[[1]][1]))
  fold_cutoff = 1
  pvalue_cutoff = 0.05
  ds$ID  <- apply(ds,1, function(x) strsplit(x[["Gene.names"]], ";")[[1]][1])
  ds$ID <- ifelse(is.na(ds$ID),ds$uniprotID,ds$ID)
  ds$ID[ds$ID=="Q9UN81"] <- "L1RE1"
  ds$ID[ds$ID=="O00370"] <- "LORF2"
  
  cold <- ifelse(grepl("ORF2p",ds$ID),"purple", ifelse(ds$ID=="L1RE1"|grepl("ORF1p",ds$ID), "red", ifelse(ds$outliers&(ds$Significant == "Yes"),"black", "grey")))
  group <- ifelse(grepl("ORF2p",ds$ID),"ORF2s", ifelse(ds$ID=="L1RE1"|grepl("ORF1p",ds$ID), "ORF1s",ifelse(ds$outliers&ds$Significant == "Yes", "Significant","Not significant/Outliers")))
  size <- ifelse((ds$Significant == "Yes")&( ds$uniprotID %in% eLife_list$uniprotID),5,1)
  ds$colour <- ifelse(grepl("ORF2p",ds$ID),"purple", ifelse(grepl("ORF1p", ds$ID)|ds$ID=="L1RE1", "red",ifelse((ds$Significant == "Yes"& ds$outliers),"black","grey")))
  lab <- ifelse((ds$Significant=="Yes"& ds$outliers)|ds$ID=="L1RE1"|grepl("ORF1p", ds$ID),ds$ID,"")
  
  cols <- c("ORF2s"="purple", "ORF1s"="red",  "Significant"="black",  "Not significant/Outliers"="grey")
  
  subds <-  subset(ds, (((ds$uniprotID %in% eLife_list$uniprotID)|ds$outliers)&Significant=="Yes")|(ID=="L1RE1")|(grepl("ORF1p",ID))|grepl("ORF2p",ID))
  #pdf(paste0("../Image/Volcano_plot/Volcano_plot_",condition,".pdf"),width = 12, height = 18)
  png(paste0("../Image/Volcano_plot/Volcano_plot_",condition,".png"),width = 1000, height = 700)
  
  p <- ggplot(ds,aes(logfold, -log10(p.adj),label=ID)) +
    geom_point(col= "orange",size=size)+
    geom_point(aes(colour=group),fill = cold, size=2) +
    geom_vline(xintercept = fold_cutoff, col = "blue")+
    geom_vline(xintercept = -fold_cutoff, col = "blue")+
    geom_hline(yintercept = -log10(pvalue_cutoff), col = "green")+
    ggtitle(condition)+theme_minimal()+
    scale_colour_manual(values=cols, aesthetics = c("fill","colour"))+
    geom_text_repel(data  = subds, 
                    colour = ifelse(grepl("ORF2p",subds$ID),"purple",ifelse(grepl("L1RE1|ORF1p",subds$ID),"red","black")), segment.size  = 0.2, 
                    segment.alpha =0.35,
                    box.padding = unit(0.45, "lines"), 
                    point.padding = unit(0.45, "lines"),
                    size=3 )
   # scale_x_continuous(limits = c(-7, 12))+
  #  scale_y_continuous(limits = c(0, 4.5))
  
  print(p)
  dev.off()
}


