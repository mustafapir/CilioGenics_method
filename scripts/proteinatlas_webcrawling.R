

library(data.table)
library(XML)
library(dplyr)
library(clipr)

proteinatlas<-fread("proteinatlas.tsv") %>%
  select(c(1,3,4:12))
proatlas1<-as.data.table(proteinatlas)

pb = txtProgressBar(min = 0, max =length(proteinatlas$Ensembl) , initial = 0)
for (i in 1:length(proteinatlas$Ensembl)){
  ensgene<-proteinatlas$Ensembl[i]
  link<-paste0("http://www.proteinatlas.org/", ensgene, ".xml")
  data <- xmlParse(link)
  xml_data<-xmlToList(data)
  cdata<-xpathSApply(data, "//text()", xmlValue)
  cdata<-gsub("\n","",cdata)
  a<-cdata[grepl("cilia", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,2:=paste(a, collapse = ",")]
    proatlas1[i,7:="YES"]
  }
  else {proatlas1[i,2:=NA]
        proatlas1[i,7:="NO"]}
  
  a<-cdata[grepl("cilium", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,3:=paste(a, collapse = ",")]
    proatlas1[i,8:="YES"]
  }
  else {proatlas1[i,3:=NA]
        proatlas1[i,8:="NO"]}
  
  a<-cdata[grepl("centrosome", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,4:=paste(a, collapse = ",")]
    proatlas1[i,9:="YES"]
  }
  else {proatlas1[i,4:=NA]
        proatlas1[i,9:="NO"]}
  
  a<-cdata[grepl("flagella", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,5:=paste(a, collapse = ",")]
    proatlas1[i,10:="YES"]
  }
  else {proatlas1[i,5:=NA]
        proatlas1[i,10:="NO"]}
  
  a<-cdata[grepl("flagellum", cdata, ignore.case = TRUE)]
  if (length(a)>0){
    proatlas1[i,1:=proteinatlas$Gene[i]]
    proatlas1[i,6:=paste(a, collapse = ",")]
    proatlas1[i,11:="YES"]
  }
  else {proatlas1[i,6:=NA]
        proatlas1[i,11:="NO"]}
  
  setTxtProgressBar(pb,i)
}
write.table(proatlas1, "proatlas1.csv", quote = FALSE, sep = "|", row.names = FALSE)


proatlas2<-fread("proatlas1.csv")
colnames(proatlas2)<-c("Gene_name", "Cilia_comments", "Cilium_comments", "Centrosome_comments", "Flagella_comments", "Flagellum_comments",
                       "Cilia", "Cilium", "Centrosome", "Flagella", "Flagellum")

proatlas3<- subset(proatlas2, Cilia == "YES" | Cilium == "YES" | Centrosome == "YES" | Flagella == "YES" | Flagellum == "YES")

write_clip(proatlas3)



m1<-read_xlsx("modified_protein_atlas_gene_list.xlsx", sheet = 2)
m2<-read_xlsx("modified_protein_atlas_gene_list.xlsx", sheet = 3)
m3<-read_xlsx("modified_protein_atlas_gene_list.xlsx", sheet = 4)
m4<-read_xlsx("modified_protein_atlas_gene_list.xlsx", sheet = 5)

m1$score<-0.25
m2$score<-0.5
m3$score<-1
m4$score<-1

protein_atlas<-rbind(m1, m2, m3, m4)
protein_atlas<-subset(protein_atlas, !is.na(Gene_name))


write.table(protein_atlas[,c(1,12)], "protein_atlas_scores.csv", row.names = FALSE, quote = FALSE)




