

`%fin%` <- function(x, table) {
  stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}



mousesynonymConverter<-function(genelist,colname){
  
  names(mouse_synonyms1)[1]<-"Gene_name_temp"
  n_occur <- data.frame(table(mouse_synonyms1$Gene_synonyms))
  bg5<-mouse_synonyms1[mouse_synonyms1$Gene_synonyms %in% n_occur$Var1[n_occur$Freq > 1],]
  
  bg5<-bg5[which(bg5$Gene_synonyms %in% genelist[[colname]][which(!(genelist[[colname]] %in% mouse_synonyms1$Gene_name))]),]
  colnames(bg5)[2]<-colname
  
  genelist<-merge(genelist, bg5, by = colname, all = TRUE, allow.cartesian = TRUE)
  for (i in 1:length(genelist[[colname]])){
    
    if (!is.na(genelist$Gene_name_temp[i])){
      genelist[[colname]][i]<-genelist$Gene_name_temp[i]
    }
  }
  pb = txtProgressBar(min = 0, max =length(genelist[[colname]]) , initial = 0)
  for (i in 1:length(genelist[[colname]])){
    if (!(genelist[[colname]][i] %fin% mouse_synonyms1$Gene_name_temp) && length(mouse_synonyms1$Gene_name_temp[which(mouse_synonyms1$Gene_synonyms %fin% genelist[[colname]][i])]) == 1){
      genelist[[colname]][i]<-mouse_synonyms1$Gene_name_temp[which(mouse_synonyms1$Gene_synonyms %fin% genelist[[colname]][i])]
    }
    setTxtProgressBar(pb,i)
  }
  return(genelist[,1:(length(genelist)-1), drop = FALSE])
}


mousegnameConverter<-function(genelist,colname){
  
  names(mouse_homology)[1]<-"Gene_name_temp"
  n_occur <- data.frame(table(mouse_homology$Gene_synonyms))
  bg5<-mouse_homology[mouse_homology$Gene_synonyms %in% n_occur$Var1[n_occur$Freq > 1],]
  
  bg5<-bg5[which(bg5$Gene_synonyms %in% genelist[[colname]][which(!(genelist[[colname]] %in% mouse_homology$Gene_name))]),]
  colnames(bg5)[2]<-"Gene_name_A"
  
  genelist<-merge(genelist, bg5, by = colname, all = TRUE, allow.cartesian = TRUE)
  for (i in 1:length(genelist[[colname]])){
    
    if (!is.na(genelist$Gene_name_temp[i])){
      genelist[[colname]][i]<-genelist$Gene_name_temp[i]
    }
  }
  
  pb = txtProgressBar(min = 0, max =length(genelist[[colname]]) , initial = 0)
  for (i in 1:length(genelist[[colname]])){
    if (!(genelist[[colname]][i] %fin% mouse_homology$Gene_name_temp) && length(mouse_homology$Gene_name_temp[which(mouse_homology$Gene_synonyms %fin% genelist[[colname]][i])]) == 1){
      genelist[[colname]][i]<-mouse_homology$Gene_name_temp[which(mouse_homology$Gene_synonyms %fin% genelist[[colname]][i])]
    }
    setTxtProgressBar(pb,i)
  }
  return(genelist[,1:(length(genelist)-1), drop = FALSE])
}






gnameConverter<-function(genelist,colname){
  
  names(gene_synonyms2)[1]<-"Gene_name_temp"
  n_occur <- data.frame(table(gene_synonyms2$Gene_synonyms))
  bg5<-gene_synonyms2[gene_synonyms2$Gene_synonyms %in% n_occur$Var1[n_occur$Freq > 1],]
  
  bg5<-bg5[which(bg5$Gene_synonyms %in% genelist[[colname]][which(!(genelist[[colname]] %in% gene_synonyms2$Gene_name))]),]
  colnames(bg5)[2]<-colname
  
  genelist<-merge(genelist, bg5, by = colname, all = TRUE, allow.cartesian = TRUE)
  for (i in 1:length(genelist[[colname]])){
    
    if (!is.na(genelist$Gene_name_temp[i])){
      genelist[[colname]][i]<-genelist$Gene_name_temp[i]
    }
  }
  
  pb = txtProgressBar(min = 0, max =length(genelist[[colname]]) , initial = 0)
  for (i in 1:length(genelist[[colname]])){
    if (!(genelist[[colname]][i] %in% gene_synonyms2$Gene_name_temp) && length(gene_synonyms2$Gene_name_temp[which(gene_synonyms2$Gene_synonyms %in% genelist[[colname]][i])]) == 1){
      genelist[[colname]][i]<-gene_synonyms2$Gene_name_temp[which(gene_synonyms2$Gene_synonyms %in% genelist[[colname]][i])]
    }
    setTxtProgressBar(pb,i)
  }
  return(genelist[,1:(length(genelist)-1), drop = FALSE])
}



apiToGene<- function(x){
  dt<-data.frame(t(x), stringsAsFactors = FALSE)
  dt<-data.frame(dt[c(-1,-2),1], stringsAsFactors = FALSE)
  dt<-str_split(dt[,1], pattern = '"symbol":"')
  dt<-data.frame(unlist(dt), stringsAsFactors = FALSE)
  dt<-data.frame(dt[seq(2, length(dt[,1])-2, 2),1], stringsAsFactors = FALSE)
  colnames(dt)[1]<-"gene_name"
  dt
}

normalization<-function(x){
  a<-(x-min(x))/(max(x)-min(x))
}


stdize = function(x, ...) {((x - min(x, ...)) / (max(x, ...) - min(x, ...)))-1}

