## Mehrnoosh Oghbaie
## 06/13/2019
## ANOVA test between cases and controls

##  1. Filter proteins with less than two peptides (I'll skip this one for now)
##  2. Perform t-test between cases and controls
##  3. Adjust pvalueswith Benjamin Hochberg correction test
##  4. Select significant proteins with (p.adj < 0.05 & log2fold > 1)
Template$set("public","experimentPreImputateSignificant", list())
Template$set("public","experimentImputateSignificant", list())


Template$set("public","anovaAnalysisPreImpute", function(){
  for(name in names(self[["PreImputation"]])){
    y <- self[["PreImputation"]][[name]]
    x <- strsplit(name,".vs.")[[1]]
    
    dt <- data.frame(matrix(NA,ncol=0,nrow=dim(y)[1]))
    dt[["ID"]] <- y[["id"]]
    dt[["uniprotID"]] <- y[["uniprotID"]]
    dt[["Gene.names"]] <- as.character(y[["Gene.names"]])
    dt[["p.value"]] <- NA
    dt[["logfold"]] <- NA
    dt[["Significant"]] <- NA
    dt[["Control_non_zero"]] <- apply(y[,colnames(y)[grepl(x[2],colnames(y))]],1, function(z) sum(!is.na(unname(unlist(z)))))
    dt[["Case_non_zero"]] <- apply(y[,colnames(y)[grepl(x[1],colnames(y))]],1, function(z) sum(!is.na(unname(unlist(z)))))
    w <- y[,-c(1:3)]
    w[is.na(w)] <- 0
    for(k in 1:dim(dt)[1]){
      print(k)
      dt[[k,"p.value"]] <- t.test(unname(unlist(w[k,grepl(x[1],colnames(w))])),unname(unlist(w[k,grepl(x[2],colnames(w))])))$p.value
      dt[[k,"logfold"]] <- mean(unname(unlist(w[k,grepl(x[1],colnames(w))])))-mean(unname(unlist(w[k,grepl(x[2],colnames(w))])))
    }
    
    dt[["p.adj"]] <- p.adjust(dt[["p.value"]], method = "BH", n = length(dt[["p.value"]]))
    dt[["Significant"]] <- ifelse(dt[["p.adj"]]<0.05&abs(dt[["logfold"]])>1,"Yes","No")
    
    ## To be changed
    self[["experimentPreImputateSignificant"]][[name]] <-dt
    
  }
})

Template$set("public","anovaAnalysis", function(){
  for(name in names(self[["experimentImputed"]])){
    y <- self[["experimentImputed"]][[name]]
    x <- strsplit(name,".vs.")[[1]]
    
    dt <- data.frame(matrix(NA,ncol=0,nrow=dim(y)[1]))
    dt[["ID"]] <- y[["id"]]
    dt[["uniprotID"]] <- y[["uniprotID"]]
    dt[["Gene.names"]] <- as.character(y[["Gene.names"]])
    dt[["p.value"]] <- NA
    dt[["logfold"]] <- NA
    dt[["Significant"]] <- NA
    dt[["Control_non_zero"]] <- apply(y[,colnames(y)[grepl(x[2],colnames(y))]],1, function(z) sum(!is.na(unname(unlist(z)))))
    dt[["Case_non_zero"]] <- apply(y[,colnames(y)[grepl(x[1],colnames(y))]],1, function(z) sum(!is.na(unname(unlist(z)))))
    w <- y[,-c(1:3)]
   # w[is.na(w)] <- 0
    for(k in 1:dim(dt)[1]){
      print(k)
      dt[[k,"p.value"]] <- t.test(unname(unlist(w[k,grepl(x[1],colnames(w))])),unname(unlist(w[k,grepl(x[2],colnames(w))])))$p.value
      dt[[k,"logfold"]] <- mean(unname(unlist(w[k,grepl(x[2],colnames(w))])))-mean(unname(unlist(w[k,grepl(x[1],colnames(w))])))
    }
    
    dt[["p.adj"]] <- p.adjust(dt[["p.value"]], method = "BH", n = length(dt[["p.value"]]))
    dt[["Significant"]] <- ifelse(dt[["p.adj"]]<0.05&abs(dt[["logfold"]])>1,"Yes","No")
    dt[!is.na(dt[["p.value"]]),"Significant_fdr_0.05"] <- qvalue(p = dt[["p.value"]][!is.na(dt[["p.value"]])], lfdr=TRUE,fdr.level = 0.05)$significant
    dt[["Significant_fdr_0.05_logfold"]] <- ifelse(dt[["Significant_fdr_0.05"]]&abs(dt[["logfold"]])>1,"Yes","No")
    dt[["lfdr"]] <- lfdr(p = dt[["p.value"]])
    
    ## To be changed
    self[["experimentImputateSignificant"]][[name]] <-dt
    
  }
})


####################################################################################################################
## Venn diagram
Template$set("public","Significant_list", NA)
Template$set("public","Significant_proteins", NA)

Template$set("public","drawVenDiagramSignificant", function(){
  Significant_proteins <- data.frame(matrix(NA, nrow=length(self[["input_merged"]]$uniprotID), ncol= 3+length( names(self[["experimentImputateSignificant"]]))))
  colnames(Significant_proteins) <- c("id","uniprotID","Gene.names", names(self[["experimentImputateSignificant"]]))
  Significant_proteins$uniprotID <- self[["input_merged"]]$uniprotID
  Significant_proteins$id <- self[["input_merged"]]$id
  Significant_proteins$Gene.names <- as.character(self[["input_merged"]]$Gene.names)
  
  for(i in names(self[["experimentImputateSignificant"]])){
    case1more <- self[["experimentPreImputateSignificant"]][[i]][["Case_non_zero"]][match(Significant_proteins$uniprotID,self[["experimentPreImputateSignificant"]][[i]][,"uniprotID"])]
    Significant_proteins[[i]] <- ifelse(case1more>1, 
                                        self[["experimentImputateSignificant"]][[i]][,"Significant"][match(Significant_proteins$uniprotID,self[["experimentImputateSignificant"]][[i]][,"uniprotID"])],
                                        "No")
    
    Significant_proteins[[i]][is.na(Significant_proteins[[i]])] <- "NO"
    Significant_proteins[[i]] <- ifelse(Significant_proteins[[i]]=="Yes",1,0)
  }
  self[["Significant_proteins"]] <- Significant_proteins
  
  Significant_list <- Significant_proteins%>%
    dplyr::filter(.[[4]]|.[[5]]|.[[6]]|.[[7]]|.[[8]]|.[[9]]==1)%>%
    dplyr::select(uniprotID, Gene.names)
  self[["Significant_list"]] <- Significant_list
  png(paste0("../Image/Volcano_plot/Significant_Venndiagram.png"),width = 800, height = 600)
  venn(Significant_proteins[,-c(1:3)],zcolor = "style")
  dev.off()
})

####################################################################################################################
## 5. Volcano plot "../Image/Volcano_plot/"
####################################################################################################################
Template$set("public","drawVolcanoPlot", function(){
  for(condition in names(PA[["experimentImputateSignificant"]])){
    draw_volcanoplot(data=PA[["experimentImputateSignificant"]],condition=condition)
  }
  png(paste0("../Image/Volcano_plot/p.value_distribution.png"),width = 1400, height = 700)
  par(mfrow=c(2,3),mar = c(2,2,2,2) +3)
  for(nam in names(PA[["experimentImputateSignificant"]])){
    p.value <- PA[["experimentImputateSignificant"]][[nam]]$p.value
    p.adj <- PA[["experimentImputateSignificant"]][[nam]]$p.adj
    hist(p.value, breaks=50, xlim=c(0,1), ylab ="p.value", xlab="", main=nam, cex.lab=1.4)
    par(new = TRUE)
    plot(density(p.adj), col="red",axes=F, xlab="", ylab="",main="", xlim=c(0,1))+
      axis(4)
    mtext("p.adj",4,2 ,col="red", cex=1)
  }
  dev.off()
})


