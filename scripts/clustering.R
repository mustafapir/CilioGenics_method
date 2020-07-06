library(readxl)
library(data.table)
library(dplyr)
library(pbapply)
library(tidyverse)
library(myTAI)
library(cluster)
library(fpc)

file<-read.csv("merged_filtered_blast1.txt", sep = "\t")
species_name<-read_excel("species_name.xlsx", col_names = FALSE)
ncbi_np_new<- read.csv("ncbi_protein_gene_31.12.19.txt", sep = "\t", stringsAsFactors = FALSE)
ncbi_npU_new<- unique(setDT(ncbi_np_new), by= "Protein.product")

orderoforg<-read_xlsx("Backup of Ciliated and Non-ciliary Organisms.xlsx")
orderoforg<-data.frame(select(orderoforg, Organism), stringsAsFactors = FALSE)
cluster_scores<-fread("cluster_scores.txt")

#ind_odd<-seq(3, 145, 2)
file1<-file

cnames<-data.frame(colnames(file1), stringsAsFactors = FALSE)
colnames(file1)<-gsub("\\.", "_", colnames(file1))
df<-data.frame(species_name$...3, stringsAsFactors = FALSE)
df1<-data.frame(species_name$...1, stringsAsFactors = FALSE)
for (i in 1:length(df$species_name....3)){
  colnames(file1)<-gsub(paste(df[i,]), paste(df1[i,]), colnames(file1))
}
colnames(file1)<-gsub("GCF_000002235_5_Spur_5_0", "Strongylocentrotus purpuratus", colnames(file1))
colnames(file1)<-gsub("GCF_000002595_1_v3_0", "Chlamydomonas reinhardtii", colnames(file1))
colnames(file1)<-gsub("GCF_000004195_4_UCB_Xtro_10_0", "Xenopus tropicalis", colnames(file1))
colnames(file1)<-gsub("GCF_000025565_1_ASM2556v1", "Enterobacter cloacae", colnames(file1))
colnames(file1)<-gsub("GCF_000090985_2_ASM9098v2", "Micromonas commoda", colnames(file1))
colnames(file1)<-gsub("GCF_008122165_1_Kamilah_GGO_v0", "Gorilla gorilla gorilla", colnames(file1))

colnames(file1)<-gsub("__", "_", colnames(file1))
colnames(file1)<-gsub(" ", "_", colnames(file1))

file2<-file1

file2$gene_name<-pblapply(file2[,1], function(x) ncbi_npU_new$Locus[match(x, ncbi_npU_new$Protein.product)])
file2<-file2[c(146, 1:145)]
file2<-file2 %>% mutate(gene_name = sapply(gene_name, toString))



odd_index2<-seq(4, 146, 2)
nscores<- select(file2, gene_name, human_np, colnames(file2[,odd_index2]))
nscores[,2]<-as.character(nscores$human_np)
nscores<-unique(setDT(nscores), by= "gene_name")

nscores2<-nscores
rownames(nscores2)<-nscores2$gene_name
nscores2[,3:74][is.na(nscores2[,3:74])] <- 0
nscores2[,3:74] <- lapply(nscores2[,3:74], factor)

colnames(nscores2)<-gsub("_e_value", "", colnames(nscores2))
orderoforg$Organism<-gsub(" ", "_", orderoforg$Organism)

colnames(nscores2)[74]<-"Homo_sapiens"
orderoforg$Organism[c(32,42,53,69)]<-NA
orderoforg<-data.frame(na.omit(orderoforg), stringsAsFactors = FALSE)

setcolorder(nscores2, orderoforg$Organism)
nscores21<-gnameConverter(nscores2, "gene_name")
nscores21<-nscores21[-12842,]
row.names(nscores21)<-nscores21$gene_name

write.table(nscores21, "nscores2.csv", quote = FALSE, sep = "\t", row.names = FALSE)


# ***********************************************************************************

nscores3<-nscores2 %>% mutate_if(is.numeric, ~1 * (. > 50))
nscores3[,1:72] <- lapply(nscores3[,1:72], factor)


species<-read_xlsx("Backup of Ciliated and Non-ciliary Organisms.xlsx") %>%
  select(Organism, Class) %>%
  separate("Class", c("Class", "Family", "class2"), sep = ";", fill = "right") %>%
  select(Organism, Family) %>%
  lapply(function(x) gsub(" ", "_", x)) %>%
  unlist() %>%
  matrix(nrow = 76, ncol = 2) %>%
  data.frame(stringsAsFactors = FALSE)

species<-species[which(species$X1 %in% colnames(nscores2[,1:72])),]
write.table(species, "species_app.txt", quote = FALSE)
anot<-data.frame(Class = species[,2], organisms = rep(c("Ciliary", "Nonciliary"), c(44,28)))
row.names(anot)<-species$X1

my_colour = list(organisms = c(Ciliary = "firebrick3", Nonciliary = "dodgerblue3"), 
                 Class = c(Animals = "firebrick3", Fungi = "dodgerblue3", Protists = "darkgrey", Plants = "chartreuse", Other = "ghostwhite", Bacteria = "gray0"))



nscores3<-nscores2
nscores3[,1:72] <- lapply(nscores3[,1:72], factor)

gower_dist<-daisy(nscores3[,1:72], metric = "gower")
aggl.clust <- hclust(gower_dist, method = "ward.D2")

tree<-cutree(aggl.clust, k=60)
aa<-cbind(nscores3[,73], tree)

write.csv()
#aa<-cbind(nscores3[,73], cutree(aggl.clust, k=60))

aa$cluster_scores<-cluster_scores$score[match(aa$tree, cluster_scores$cluster_number)]
aa$norm_cluster_scores<-normalization(aa$cluster_scores)
colnames(aa)[1]<-"Gene_name"

write.table()
aa<-read.table("aa.txt", stringsAsFactors = FALSE)

gnameConverter(aa$gene_name)

colnames(aa)<-c("Gene_name", "cluster_number")
clstr_with_scores<-merge(aa, cluster_scores, by = "cluster_number")
gnameConverter(protAtlas$gene_symbol)

scores<-nscores[1:60,1:2]
for (i in 1:60){
  scores$cluster[i]<-i
  scores$ciliary[i]<-length(na.omit(match(ciliaryGenes1$Gene.Name, aa$gene_name[which(aa$tree==i)])))
  scores$all[i]<-length(which(aa$tree==i))
  scores$nonciliary[i]<-length(na.omit(match(not_ciliary$`Gene Name`, aa$gene_name[which(aa$tree==i)])))
  scores$score[i]<-(length(na.omit(match(ciliaryGenes1$Gene.Name, aa$gene_name[which(aa$tree==i)])))-
                      length(na.omit(match(not_ciliary$`Gene Name`, aa$gene_name[which(aa$tree==i)]))))/length(which(aa$tree==i))
  
}

cluster_scores

protAtlas<-read_xlsx("protein_atlas2.xlsx")
proAtlasScores<-ord_merged_scores3[which(ord_merged_scores3$gene_name %in% protAtlas$gene_symbol),]



write.csv(known, file = "known.csv", quote = FALSE)
write.csv(known_ciliary, file = "known_ciliary.csv", quote = FALSE)
write.csv(known_not_ciliary, file = "known_not_ciliary.csv", quote = FALSE)
write.csv(ord_merged_scores3, file = "ord_merged_scores3.csv", quote = FALSE)




