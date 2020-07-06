

library(readxl)
library(dplyr)
library(data.table)
library(pbapply)
library(stringr)
library(tidyr)


source("functions.R")

# Load the necessary files

gene_synonyms<-fread("Homo_sapiens.gene_info")
gene_synonyms2<-select(gene_synonyms, Symbol, Synonyms)
s<-strsplit(as.character(gene_synonyms2$Synonyms), '\\|')
gene_synonyms2<-data.frame(Gene_name = rep(gene_synonyms2$Symbol, sapply(s, FUN=length)), Gene_synonyms=unlist(s), stringsAsFactors = FALSE)

mouse_synonyms<-fread("Mus_musculus.gene_info")
mouse_synonyms<-select(mouse_synonyms, Symbol, Synonyms)
s<-strsplit(as.character(mouse_synonyms$Synonyms), '\\|')
mouse_synonyms<-data.frame(Gene_name = rep(mouse_synonyms$Symbol, sapply(s, FUN=length)), Gene_synonyms=unlist(s), stringsAsFactors = FALSE)
mouse_synonyms1<-mouse_synonyms[which(!(mouse_synonyms$Gene_synonyms %in% gene_synonyms2$Gene_synonyms)),]
mouse_synonyms1<-mouse_synonyms1[which(!(mouse_synonyms1$Gene_synonyms %in% gene_synonyms2$Gene_name)),]

mouse_homology<- fread("mouse_human_homology_all.csv")
mouse_homology<-mouse_homology[,-1]
mouse_homology<-mouse_homology[,-2]
mouse_homology<-unique(mouse_homology)
names(mouse_homology)<-c("Gene_name", "Gene_synonyms")
mouse_homology1<-mouse_homology[which(!(mouse_homology$Gene_synonyms %in% gene_synonyms2$Gene_name)),]
mouse_homology1<-mouse_homology1[which(!(mouse_homology1$Gene_synonyms %in% gene_synonyms2$Gene_synonyms)),]
mouse_synonyms1<-mouse_synonyms1[which(!(mouse_synonyms1$Gene_synonyms %in% mouse_homology1$Gene_synonyms)),]

mouse_homology1<-mousesynonymConverter(mouse_homology, "Gene_synonyms")
mouse_homology<-select(mouse_homology1, Gene_name, Gene_synonyms)

cluster_scores<-fread("cluster_scores.txt")

ncbi_np_new<- read.csv("ncbi_protein_gene_31.12.19.txt", sep = "\t", stringsAsFactors = FALSE)
not_ciliary<-read_xls("Nevers_2017_NegativeGenesInsightsCiiaryGenes_SuppTable3 .xls")
not_ciliary<-not_ciliary[,1]
not_ciliary$scores<-1
colnames(not_ciliary)[1]<-"Gene.Name"
not_ciliary<-gnameConverter(not_ciliary, "Gene.Name")
negative<-not_ciliary
not_ciliary<-data.frame(not_ciliary, stringsAsFactors = FALSE)
negative<-data.frame(negative, stringsAsFactors = FALSE)


ciliaryGenes<-read.csv("OurOwn_NewGoldStandartList_2019.tsv", sep = "\t", stringsAsFactors = FALSE)
ciliaryGenes1<-select(ciliaryGenes, Gene.Name)
ciliaryGenes1<-unique(ciliaryGenes1)
ciliaryGenes1[433,]<-"RABL2A"
ciliaryGenes1$scores<-1

ciliaryGenes1<-gnameConverter(ciliaryGenes1, "Gene.Name")


reyfman<-read.table("Reyfman_Table1.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
lung <-read_xlsx("Travaglini_KJ_TablesS4_Human_lung_Single_cell_RNA_seq.xlsx", sheet = 3, skip = 1)
lung1<-read_xlsx("Travaglini_KJ_TablesS4_Human_lung_Single_cell_RNA_seq.xlsx", sheet = 4, skip = 1)
lung2 <-read_xlsx("Travaglini_KJ_TablesS4_Human_lung_Single_cell_RNA_seq.xlsx", sheet = 5, skip = 1)

mouse<-read.table("cilia_mouse_human_homology.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
mouse<-data.frame(select(mouse, human_name), stringsAsFactors = FALSE)
colnames(mouse)[1]<-"Gene_name"

mouse<-gnameConverter(mouse, "Gene_name")

c_elegans<-read.table("C_elegans_single_cell.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
c_elegans1<-data.frame(select(c_elegans, symbol, Ciliated_sensory_neurons), stringsAsFactors = FALSE)
homology_celegans<-read_excel("celegans_homology.xlsx") %>%
  unique()

colnames(homology_celegans)[1]<-"gene_id"
celegans_hom<-fread("celegans_homology.txt")
celegans_hom<-celegans_hom[-which(celegans_hom$`Human gene name`=="")]

biogrid<-fread("BIOGRID-ORGANISM-Homo_sapiens-3.5.180.tab2.txt", sep = "\t", stringsAsFactors = FALSE)
intact<-fread("intact.txt", sep = "\t", stringsAsFactors = FALSE, quote = "")
gname_unipro<-fread("uniprot_gname.tab", sep = "\t", stringsAsFactors = FALSE)

wormbase<-read.delim("c_elegans.PRJNA13758.WS274.interactions.txt", sep = "\t", colClasses = c(rep("character", 11), rep("NULL", 78)), header = FALSE,
                     col.names = paste0("V",seq_len(89)), fill = TRUE, na.strings = c(""," ","NA"))
wormbase<-wormbase[-1,]
wormbase<-wormbase[-1,]
wormbase<-wormbase[-1,]
colnames(wormbase)<-wormbase[1,]
wormbase<-wormbase[-1,]


LM1_RFX1<-fread("LM1_RFX1") %>% select(gene_name) %>% unique()
RFX1<-fread("RFX1_M00280") %>% select(gene_name) %>% unique()
RFX1_2<-fread("RFX1_M00281") %>% select(gene_name) %>% unique()
foxj1<-fread("foxj1.txt", sep = ",", header = FALSE)
rfx1_transfac<-fread("RFX1_TRANSFAC.txt", sep = ",", header = FALSE)

foxj1<-apiToGene(foxj1)
rfx1_transfac<-apiToGene(rfx1_transfac)

all_pub2<-fread("all_publications.txt", stringsAsFactors = FALSE)

# Single cell

lung<- data.frame(lung$Gene, stringsAsFactors = FALSE)
lung1<- data.frame(lung1$Gene, stringsAsFactors = FALSE)
lung2<- data.frame(lung2$Gene, stringsAsFactors = FALSE)

colnames(lung)[1]<-"Gene_name"
colnames(lung1)[1]<-"Gene_name"
colnames(lung2)[1]<-"Gene_name"


lung_cilia<-rbind(lung, lung1)
lung_cilia<-rbind(lung_cilia, lung2)
lung_cilia1<-data.frame(unique(lung_cilia$Gene_name), stringsAsFactors = FALSE)
colnames(lung_cilia1)<-"Gene_name"

pb = txtProgressBar(min = 1, max = 88, initial = 0)
for (i in c(1:2,6:88)){
  file1<-read_xlsx("Travaglini_KJ_TablesS4_Human_lung_Single_cell_RNA_seq.xlsx", sheet = i, skip = 1)
  colnames(file1)[1]<-"Gene_name"
  file1<-data.frame(file1$Gene_name, stringsAsFactors = FALSE)
  colnames(file1)[1]<-"Gene_name"
  lung_cilia1<-data.frame(lung_cilia1$Gene_name[which(!(lung_cilia1$Gene_name %in% file1$Gene_name))], stringsAsFactors = FALSE)
  colnames(lung_cilia1)<-"Gene_name"
  print(c(length(lung_cilia1$Gene_name), i))
  setTxtProgressBar(pb,i)
}

lung_cilia1<-gnameConverter(lung_cilia1, "Gene_name")


reyfman1<-select(reyfman, Ciliated.Cells)
colnames(reyfman1)<-"Gene_name"
for (i in 1:13){
  reyfman1<-data.frame(reyfman1$Gene_name[which(!(reyfman1$Gene_name %in% reyfman[[1]]))], stringsAsFactors = FALSE)
  colnames(reyfman1)<-"Gene_name"
}

reyfman1<-gnameConverter(reyfman1, "Gene_name")

human<-rbind(reyfman1, lung_cilia1)
human$score<-1
human$score[duplicated(human$Gene_name)]<-2
human1<-human[order(-human$score),]
human1<-unique(setDT(human1), by= "Gene_name")


# C. elegans single cell

c_ele<- c_elegans
pb = txtProgressBar(min = 0, max = length(c_elegans$symbol), initial = 0)
for (i in 1:length(c_elegans$symbol)){
  if (length(which(c_elegans[i, c(3:5, 7:29)]>15))<4 && c_elegans[i,6]>15){
    
    if (all(c_elegans[i,c(3:5, 7:19, 21:29)] < c_elegans[i,6])){
      c_ele[i,]<-c_elegans[i,]
    }
    else{
      c_ele[i,]<-NA
    }
  }
  else{
    c_ele[i,]<-NA
  }
  setTxtProgressBar(pb,i)
}
c_ele<-na.omit(c_ele)

new_cele<-merge(c_ele, homology_celegans, by="gene_id")
new_cele<-na.omit(new_cele)
colnames(new_cele)[30]<-"Gene_name"
new_cele2<-select(new_cele, Gene_name)
new_cele2<-gnameConverter(new_cele2, "Gene_name")


# Protein - Protein interactions


# Intact #
intact0<-intact[-which(!grepl("taxid:9606", intact$`Taxid interactor A`))]
intact0<-intact0[-which(!grepl("taxid:9606", intact0$`Taxid interactor B`))]

intact1<-intact0[-which(grepl("EBI", intact0$`#ID(s) interactor A`))]
intact1<-intact1[-which(grepl("EBI", intact1$`ID(s) interactor B`))]

intact2<-intact1[-which(!grepl("unipro", intact1$`ID(s) interactor B`))]
intact2<-intact2[-which(!grepl("unipro", intact2$`#ID(s) interactor A`))]

intact2$`#ID(s) interactor A`<- gsub("uniprotkb:","\\2",intact2$`#ID(s) interactor A`)
intact2$`ID(s) interactor B`<- gsub("uniprotkb:","\\2",intact2$`ID(s) interactor B`)
intact2$`#ID(s) interactor A`<- gsub("-.", "\\1", intact2$`#ID(s) interactor A`)
intact2$`ID(s) interactor B`<- gsub("-.","\\1",intact2$`ID(s) interactor B`)
intact2_ids<-intact2[,1:2]


colnames(gname_unipro)[1:2]<-c(colnames(intact2_ids)[1], "Gene_name_A")
intact_names<-merge(intact2_ids, gname_unipro, by="#ID(s) interactor A")
colnames(gname_unipro)[1:2]<-c(colnames(intact2_ids)[2], "Gene_name_B")
intact_names<-merge(intact_names, gname_unipro, by="ID(s) interactor B")
intact_names<-unique(intact_names)
intact_gnames<-intact_names[,3:4]

# Biogrid #
names_biogrid<-biogrid[,8:9]
biogrid_gnames<-unique(names_biogrid)
colnames(biogrid_gnames)<-c("Gene_name_A", "Gene_name_B")

biogrid_gnames<-data.frame(biogrid_gnames, stringsAsFactors = FALSE)

biogrid_gnames<-mousesynonymConverter(biogrid_gnames,"Gene_name_A")
biogrid_gnames<-mousesynonymConverter(biogrid_gnames,"Gene_name_B")
biogrid_gnames<-mousegnameConverter(biogrid_gnames,"Gene_name_A")
biogrid_gnames<-mousegnameConverter(biogrid_gnames,"Gene_name_B")
biogrid_gnames<-gnameConverter(biogrid_gnames, "Gene_name_A")
biogrid_gnames<-gnameConverter(biogrid_gnames, "Gene_name_B")


intact_gnames<-gnameConverter(intact_gnames, "Gene_name_A")
intact_gnames<-gnameConverter(intact_gnames, "Gene_name_B")

biogrid_gnames$commons<-1
biogrid_gnames$commons[biogrid_gnames$Gene_name_A %fin% intact_gnames$Gene_name_A & biogrid_gnames$Gene_name_B %fin% intact_gnames$Gene_name_B]<-2
intact_gnames$commons<-1
intact_gnames$commons[intact_gnames$Gene_name_A %fin% biogrid_gnames$Gene_name_A & intact_gnames$Gene_name_B %fin% biogrid_gnames$Gene_name_B]<-2


# Scoring single cell

x<-reyfman1
y<-lung_cilia1

a<-length(na.omit(match(ciliaryGenes1$Gene.Name, x$Gene_name)))/length(x$Gene_name)
a1<-length(na.omit(match(ciliaryGenes1$Gene.Name, y$Gene_name)))/length(y$Gene_name)
b<-length(na.omit(match(ciliaryGenes1$Gene.Name, lung_cilia1$Gene_name)))/length(lung_cilia1$Gene_name)
c<-length(na.omit(match(ciliaryGenes1$Gene.Name, mouse$Gene_name)))/length(mouse$Gene_name)
d<-length(na.omit(match(ciliaryGenes1$Gene.Name, new_cele2$Gene_name)))/length(new_cele2$Gene_name)

a2<-length(na.omit(match(negative$Gene.Name, x$Gene_name)))/length(x$Gene_name)
b2<-length(na.omit(match(negative$Gene.Name, y$Gene_name)))/length(y$Gene_name)
c2<-length(na.omit(match(negative$Gene.Name, mouse$Gene_name)))/length(mouse$Gene_name)
d2<-length(na.omit(match(negative$Gene.Name, new_cele2$Gene_name)))/length(new_cele2$Gene_name)


singleCellscore<-data.frame(c(reyfman1, a-a2), stringsAsFactors = FALSE)
lungscore<-data.frame(c(lung, b-b2), stringsAsFactors = FALSE)
mousescore<-data.frame(c(mouse, c-c2), stringsAsFactors = FALSE)
celescore<-data.frame(c(new_cele2, d-d2), stringsAsFactors = FALSE)
singleCellscore<-merge(merge(merge(singleCellscore, lungscore, by="Gene_name", all = TRUE), mousescore, by="Gene_name", all = TRUE),
                       celescore, by="Gene_name", all = TRUE)
colnames(singleCellscore)<-c("Gene_name","Reyfman_score", "Lung_score", "Mouse_score", "C.elegans_score")
singleCellscore[is.na(singleCellscore)]<-0
singleCellscore$Single_cell_score<-rowSums(singleCellscore[,-1])

singleCellscore$normalized_values<-(singleCellscore$Single_cell_score- min(singleCellscore$Single_cell_score))/(max(singleCellscore$Single_cell_score)-min(singleCellscore$Single_cell_score))
singleCellscore_u<-unique(singleCellscore)


# Biogrid scoring #
biogrid_gnames$score<-0
biogrid_gnames$score<-pblapply(biogrid_gnames$Gene_name_B, function(x) ciliaryGenes1$scores[match(x, ciliaryGenes1$Gene.Name)])
biogrid_gnames$score2<-0
biogrid_gnames$score2<-pblapply(biogrid_gnames$Gene_name_A, function(x) ciliaryGenes1$scores[match(x, ciliaryGenes1$Gene.Name)])
biogrid_gnames$notscore<-0
biogrid_gnames$notscore<-pblapply(biogrid_gnames$Gene_name_B, function(x) not_ciliary$scores[match(x, not_ciliary$Gene.Name)])
biogrid_gnames$notscore2<-0
biogrid_gnames$notscore2<-pblapply(biogrid_gnames$Gene_name_A, function(x) not_ciliary$scores[match(x, not_ciliary$Gene.Name)])


biogrid_gnames$score[which(is.na(biogrid_gnames$score))]<-0
biogrid_gnames$score2[which(is.na(biogrid_gnames$score2))]<-0
biogrid_gnames$notscore[which(is.na(biogrid_gnames$notscore))]<-0
biogrid_gnames$notscore2[which(is.na(biogrid_gnames$notscore2))]<-0
biogrid_gnames$score<-as.numeric(biogrid_gnames$score)
biogrid_gnames$score2<-as.numeric(biogrid_gnames$score2)
biogrid_gnames$notscore<-as.numeric(biogrid_gnames$notscore)
biogrid_gnames$notscore2<-as.numeric(biogrid_gnames$notscore2)

biogrid_gnames$score<-biogrid_gnames$score*biogrid_gnames$commons
biogrid_gnames$score2<-biogrid_gnames$score2*biogrid_gnames$commons
biogrid_gnames$notscore<-biogrid_gnames$notscore*biogrid_gnames$commons
biogrid_gnames$notscore2<-biogrid_gnames$notscore2*biogrid_gnames$commons


# Intact scoring #
intact_gnames$score<-0
intact_gnames$score<-pblapply(intact_gnames$Gene_name_B, function(x) ciliaryGenes1$scores[match(x, ciliaryGenes1$Gene.Name)])
intact_gnames$score2<-0
intact_gnames$score2<-pblapply(intact_gnames$Gene_name_A, function(x) ciliaryGenes1$scores[match(x, ciliaryGenes1$Gene.Name)])
intact_gnames$notscore<-0
intact_gnames$notscore<-pblapply(intact_gnames$Gene_name_B, function(x) not_ciliary$scores[match(x, not_ciliary$Gene.Name)])
intact_gnames$notscore2<-0
intact_gnames$notscore2<-pblapply(intact_gnames$Gene_name_A, function(x) not_ciliary$scores[match(x, not_ciliary$Gene.Name)])

intact_gnames$score[which(is.na(intact_gnames$score))]<-0
intact_gnames$score2[which(is.na(intact_gnames$score2))]<-0
intact_gnames$notscore[which(is.na(intact_gnames$notscore))]<-0
intact_gnames$notscore2[which(is.na(intact_gnames$notscore2))]<-0
intact_gnames$score<-as.numeric(intact_gnames$score)
intact_gnames$score2<-as.numeric(intact_gnames$score2)
intact_gnames$notscore<-as.numeric(intact_gnames$notscore)
intact_gnames$notscore2<-as.numeric(intact_gnames$notscore2)

intact_gnames$score<-intact_gnames$score*intact_gnames$commons
intact_gnames$score2<-intact_gnames$score2*intact_gnames$commons
intact_gnames$notscore<-intact_gnames$notscore*intact_gnames$commons
intact_gnames$notscore2<-intact_gnames$notscore2*intact_gnames$commons


# Biogrid #

genlist_biogrid<-data.frame(Gene_name_A = unique(unlist(biogrid_gnames[,1:2])), stringsAsFactors = FALSE)


pb = txtProgressBar(min = 0, max = length(genlist_biogrid$Gene_name_A), initial = 0)
for (i in 1:length(genlist_biogrid$Gene_name_A)){
  x<-biogrid_gnames$score[biogrid_gnames$Gene_name_A %fin% genlist_biogrid$Gene_name_A[i]]
  x1<-biogrid_gnames$score2[biogrid_gnames$Gene_name_B %fin% genlist_biogrid$Gene_name_A[i]]
  y<-biogrid_gnames$notscore[biogrid_gnames$Gene_name_A %fin% genlist_biogrid$Gene_name_A[i]]
  y1<-biogrid_gnames$notscore2[biogrid_gnames$Gene_name_B %fin% genlist_biogrid$Gene_name_A[i]]
  
  
  genlist_biogrid$score[i]<-(length(which(x>0))+length(which(x1>0)))/(length(x)+length(x1))
  genlist_biogrid$notscore[i]<-(length(which(y>0))+length(which(y1>0)))/(length(x)+length(x1))
  genlist_biogrid$truescore[i]<-(sum(x)+sum(x1))/(length(x)+length(x1))
  if (length(x)<6){
    genlist_biogrid$scorepen2[i]<-((sum(x))-sum(y))/length(x)
  }
  else if (length(x)>5 && length(x)<11){
    genlist_biogrid$scorepen2[i]<-3*((sum(x))-sum(y))/length(x)
  }
  else if (length(x)>10 && length(x)<21){
    genlist_biogrid$scorepen2[i]<-4*((sum(x))-sum(y))/length(x)
  }
  else if (length(x)>20 && length(x)<31){
    genlist_biogrid$scorepen2[i]<-5*((sum(x))-sum(y))/length(x)
  }
  else if (length(x)>30 && length(x)<51){
    genlist_biogrid$scorepen2[i]<-6*((sum(x))-sum(y))/length(x)
  }
  else{genlist_biogrid$scorepen2[i]<-7*((sum(x))-sum(y))/length(x)}
  
  a<-length(x)+length(x1)
  
  if (a<6){
    genlist_biogrid$scorepen3[i]<-(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else if (a>5 && a<11){
    genlist_biogrid$scorepen3[i]<-3*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else if (a>10 && a<21){
    genlist_biogrid$scorepen3[i]<-4*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else if (a>20 && a<31){
    genlist_biogrid$scorepen3[i]<-5*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else if (a>30 && a<51){
    genlist_biogrid$scorepen3[i]<-6*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else{genlist_biogrid$scorepen3[i]<-7*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))}
  
  
  genlist_biogrid$newscore[i]<-(length(which(x > 0))*(((sum(x)))-sum(y)))/length(x)
  genlist_biogrid$length[i]<-length(x)
  genlist_biogrid$length1[i]<-length(x1)
  genlist_biogrid$truelength[i]<-sum(x)
  genlist_biogrid$truelength1[i]<-sum(x1)
  setTxtProgressBar(pb,i)
}


# Intact #

genlist<-data.frame(Gene_name_A = unique(unlist(intact_gnames[,1:2])), stringsAsFactors = FALSE)


pb = txtProgressBar(min = 0, max = length(genlist$Gene_name_A), initial = 0)
for (i in 1:length(genlist$Gene_name_A)){
  x<-intact_gnames$score[intact_gnames$Gene_name_A %in% genlist$Gene_name_A[i]]
  x1<-intact_gnames$score2[intact_gnames$Gene_name_B %in% genlist$Gene_name_A[i]]
  y<-intact_gnames$notscore[intact_gnames$Gene_name_A %in% genlist$Gene_name_A[i]]
  y2<-intact_gnames$notscore2[intact_gnames$Gene_name_B %in% genlist$Gene_name_A[i]]
  
  genlist$score[i]<-length(which(x > 0))/length(x)
  genlist$notscore[i]<-length(which(y > 0))/length(x)
  genlist$truescore[i]<-sum(x)/length(x)
  if (length(x)<6){
    genlist$scorepen2[i]<-((sum(x))-sum(y))/length(x)
  }
  else if (length(x)>5 && length(x)<11){
    genlist$scorepen2[i]<-3*((sum(x))-sum(y))/length(x)
  }
  else if (length(x)>10 && length(x)<21){
    genlist$scorepen2[i]<-4*((sum(x))-sum(y))/length(x)
  }
  else if (length(x)>20 && length(x)<31){
    genlist$scorepen2[i]<-5*((sum(x))-sum(y))/length(x)
  }
  else if (length(x)>30 && length(x)<51){
    genlist$scorepen2[i]<-6*((sum(x))-sum(y))/length(x)
  }
  else{genlist$scorepen2[i]<-7*((sum(x))-sum(y))/length(x)}
  
  
  if ((length(x)+length(x1))<6){
    genlist$scorepen3[i]<-(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else if ((length(x)+length(x1))>5 && (length(x)+length(x1))<11){
    genlist$scorepen3[i]<-3*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else if ((length(x)+length(x1)) && (length(x)+length(x1))<21){
    genlist$scorepen3[i]<-4*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else if ((length(x)+length(x1))>20 && (length(x)+length(x1))<31){
    genlist$scorepen3[i]<-5*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else if ((length(x)+length(x1)) && (length(x)+length(x1))<51){
    genlist$scorepen3[i]<-6*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))
  }
  else{genlist$scorepen3[i]<-7*(sum(x)+sum(x1)-sum(y)-sum(y1))/(length(x)+length(x1))}
  
  
  genlist$newscore[i]<-(length(which(x > 0))*(((sum(x)))-sum(y)))/length(x)
  genlist$length[i]<-length(x)
  genlist$length1[i]<-length(x)+length(x1)
  genlist$truelength[i]<-sum(x)
  genlist$truelength1[i]<-sum(x1)
  setTxtProgressBar(pb,i)
}

# Wormbase

wormbaseP<-wormbase[which(wormbase$Interaction_subtype =="ProteinProtein"),]
colnames(wormbaseP)[10]<-"Common-name-2"
wormbaseP<-select(wormbaseP, Interactor1, Interactor2)


wormbaseP1<-wormbaseP
wormbaseP1[]<-data.frame(pblapply(wormbaseP[,1:2], function(x) homology_celegans$Human_Symbol[match(x, homology_celegans$gene_id)]), stringsAsFactors = FALSE)
wormbaseP1<-na.omit(wormbaseP1)

wormbase_genelist<-wormbaseP1[,1]
wormbase_genelist<-data.frame(unique(wormbase_genelist), stringsAsFactors = FALSE)

wormbase_genelist2<-wormbaseP1[,2]
wormbase_genelist2<-data.frame(unique(wormbase_genelist2), stringsAsFactors = FALSE)
colnames(wormbase_genelist)<-"Gene_name"
colnames(wormbase_genelist2)<-"Gene_name"
genlist_wormbase<-rbind(wormbase_genelist, wormbase_genelist2)
genlist_wormbase<-data.frame(unique(genlist_wormbase), stringsAsFactors = FALSE)

genlist_wormbase<-gnameConverter(genlist_wormbase, "Gene_name")

wormbaseP1<-gnameConverter(wormbaseP1, "Interactor2")
wormbaseP1<-gnameConverter(wormbaseP1, "Interactor1")


wormbaseP1$score<-pblapply(wormbaseP1$Interactor2, function(x) ciliaryGenes1$scores[match(x, ciliaryGenes1$Gene.Name)])
wormbaseP1$score2<-pblapply(wormbaseP1$Interactor1, function(x) ciliaryGenes1$scores[match(x, ciliaryGenes1$Gene.Name)])
wormbaseP1$notscore<-pblapply(wormbaseP1$Interactor2, function(x) not_ciliary$scores[match(x, not_ciliary$Gene.Name)])
wormbaseP1$notscore1<-pblapply(wormbaseP1$Interactor1, function(x) not_ciliary$scores[match(x, not_ciliary$Gene.Name)])
wormbaseP1<-data.frame(matrix(unlist(wormbaseP1), nrow=length(wormbaseP1$score), ncol=length(wormbaseP1)),stringsAsFactors=FALSE)
wormbaseP1<-unique(wormbaseP1)

wormbaseP1$X3[which(is.na(wormbaseP1$X3))]<-0
wormbaseP1$X4[which(is.na(wormbaseP1$X4))]<-0
wormbaseP1$X5[which(is.na(wormbaseP1$X5))]<-0
wormbaseP1$X6[which(is.na(wormbaseP1$X6))]<-0
wormbaseP1$X3<-as.numeric(wormbaseP1$X3)
wormbaseP1$X4<-as.numeric(wormbaseP1$X4)
wormbaseP1$X5<-as.numeric(wormbaseP1$X5)
wormbaseP1$X6<-as.numeric(wormbaseP1$X6)

a<-wormbaseP1[,c(1,3)]
colnames(a)<-c("Gene_name", "score1")
b<-wormbaseP1[,c(2,4)]
colnames(b)<-c("Gene_name", "score2")
genlist_wormbase2<-merge(genlist_wormbase, a, by= "Gene_name", all = TRUE)
genlist_wormbase2<-unique(setDT(genlist_wormbase2), by="Gene_name")
genlist_wormbase2<-merge(genlist_wormbase2, b, by= "Gene_name", all = TRUE)
genlist_wormbase2<-unique(setDT(genlist_wormbase2), by="Gene_name")
genlist_wormbase2$score1[which(is.na(genlist_wormbase2$score1))]<-0
genlist_wormbase2$score2[which(is.na(genlist_wormbase2$score2))]<-0
genlist_wormbase2$total_score<-rowSums(genlist_wormbase2[,2:3])

genlist_wormbase2$total_score[which(genlist_wormbase2$total_score==2)]<-1




pb = txtProgressBar(min = 0, max = length(genlist_wormbase2$Gene_name), initial = 0)
for (i in 1:length(genlist_wormbase2$Gene_name)){
  x<-wormbaseP1$X3[wormbaseP1$X1 %fin% genlist_wormbase2$Gene_name[i]]
  y<-wormbaseP1$X5[wormbaseP1$X1 %fin% genlist_wormbase2$Gene_name[i]]
  x1<-wormbaseP1$X4[wormbaseP1$X2 %fin% genlist_wormbase2$Gene_name[i]]
  y1<-wormbaseP1$X6[wormbaseP1$X2 %fin% genlist_wormbase2$Gene_name[i]]
  if (length(x)==0){x<-0}
  if (length(y)==0){y<-0}
  if (length(x1)==0){x1<-0}
  if (length(y1)==0){y2<-0}
  
  genlist_wormbase2$truescore[i]<-(length(which(x>0))+length(which(x1>0))-length(which(y>0))-length(which(y1>0)))/(length(x)+length(x1))
  genlist_wormbase2$scorelength[i]<-length(which(x>0))+length(which(x1>0))
  genlist_wormbase2$notscorelength[i]<-length(which(y>0))+length(which(y1>0))
  genlist_wormbase2$totallength[i]<-length(x)+length(x1)
  setTxtProgressBar(pb,i)
}



# Genetic interactions

wormbaseG<-wormbase[which(wormbase$Interaction_type=="Genetic"),]
wormbaseR<-wormbase[which(wormbase$Interaction_type=="Regulatory"),]
wormbase2<-rbind(wormbaseG, wormbaseR)
colnames(wormbase2)[10]<-"Common-name-2"
wormbase2<-select(wormbase2, Interactor1, Interactor2)


# WB id to Gene name conversion
wormbase3<-wormbase2
wormbase3[]<-data.frame(pblapply(wormbase2[,1:2], function(x) homology_celegans$Human_Symbol[match(x, homology_celegans$gene_id)]), stringsAsFactors = FALSE)
wormbase3<-na.omit(wormbase3)

wormbase_glist<-wormbase3[,1]
wormbase_glist<-data.frame(unique(wormbase_glist), stringsAsFactors = FALSE)

wormbase_glist2<-wormbase3[,2]
wormbase_glist2<-data.frame(unique(wormbase_glist2), stringsAsFactors = FALSE)
colnames(wormbase_glist)<-"Gene_name"
colnames(wormbase_glist2)<-"Gene_name"
glist_wormbase<-rbind(wormbase_glist, wormbase_glist2)
glist_wormbase<-data.frame(unique(glist_wormbase), stringsAsFactors = FALSE)

wormbase3<-gnameConverter(wormbase3, "Interactor2")
wormbase3<-gnameConverter(wormbase3, "Interactor1")

glist_wormbase<-gnameConverter(glist_wormbase, "Gene_name")


wormbase3$score<-pblapply(wormbase3$Interactor2, function(x) ciliaryGenes1$scores[match(x, ciliaryGenes1$Gene.Name)])
wormbase3$score2<-pblapply(wormbase3$Interactor1, function(x) ciliaryGenes1$scores[match(x, ciliaryGenes1$Gene.Name)])
wormbase3$notscore<-pblapply(wormbase3$Interactor2, function(x) not_ciliary$scores[match(x, not_ciliary$Gene.Name)])
wormbase3$notscore1<-pblapply(wormbase3$Interactor1, function(x) not_ciliary$scores[match(x, not_ciliary$Gene.Name)])
wormbase3<-data.frame(matrix(unlist(wormbase3), nrow=length(wormbase3$score), ncol=length(wormbase3)),stringsAsFactors=FALSE)
wormbase3<-unique(wormbase3)

wormbase3$X3[which(is.na(wormbase3$X3))]<-0
wormbase3$X4[which(is.na(wormbase3$X4))]<-0
wormbase3$X5[which(is.na(wormbase3$X5))]<-0
wormbase3$X6[which(is.na(wormbase3$X6))]<-0
wormbase3$X3<-as.numeric(wormbase3$X3)
wormbase3$X4<-as.numeric(wormbase3$X4)
wormbase3$X5<-as.numeric(wormbase3$X5)
wormbase3$X6<-as.numeric(wormbase3$X6)

a<-wormbase3[,c(1,3)]
colnames(a)<-c("Gene_name", "score1")
b<-wormbase3[,c(2,4)]
colnames(b)<-c("Gene_name", "score2")
glist_wormbase2<-merge(glist_wormbase, a, by= "Gene_name", all = TRUE)
glist_wormbase2<-unique(setDT(glist_wormbase2), by="Gene_name")
glist_wormbase2<-merge(glist_wormbase2, b, by= "Gene_name", all = TRUE)
glist_wormbase2<-unique(setDT(glist_wormbase2), by="Gene_name")
glist_wormbase2$score1[which(is.na(glist_wormbase2$score1))]<-0
glist_wormbase2$score2[which(is.na(glist_wormbase2$score2))]<-0
glist_wormbase2$total_score<-rowSums(glist_wormbase2[,2:3])

glist_wormbase2$total_score[which(glist_wormbase2$total_score==2)]<-1



pb = txtProgressBar(min = 0, max = length(glist_wormbase2$Gene_name), initial = 0)
for (i in 1:length(glist_wormbase2$Gene_name)){
  x<-wormbase3$X3[wormbase3$X1 %fin% glist_wormbase2$Gene_name[i]]
  y<-wormbase3$X5[wormbase3$X1 %fin% glist_wormbase2$Gene_name[i]]
  x1<-wormbase3$X4[wormbase3$X2 %fin% glist_wormbase2$Gene_name[i]]
  y1<-wormbase3$X6[wormbase3$X2 %fin% glist_wormbase2$Gene_name[i]]
  if (length(x)==0){x<-0}
  if (length(y)==0){y<-0}
  if (length(x1)==0){x1<-0}
  if (length(y1)==0){y2<-0}
  
  glist_wormbase2$truescore[i]<-(length(which(x>0))+length(which(x1>0))-length(which(y>0))-length(which(y1>0)))/(length(x)+length(x1))
  glist_wormbase2$scorelength[i]<-length(which(x>0))+length(which(x1>0))
  glist_wormbase2$notscorelength[i]<-length(which(y>0))+length(which(y1>0))
  glist_wormbase2$totallength[i]<-length(x)+length(x1)
  setTxtProgressBar(pb,i)
}


# Prot int. score combine

score_intact<-select(genlist, Gene_name_A, scorepen3)
score_biogrid<-select(genlist_biogrid, Gene_name_A, scorepen3)
score_wb_pro<-select(genlist_wormbase2, Gene_name, truescore)
colnames(score_wb_pro)<-c("Gene_name_A", "scorepen3")
pro_score_m<-merge(score_intact, score_biogrid, by = "Gene_name_A", all = TRUE)
pro_score_m<-merge(pro_score_m, score_wb_pro, by = "Gene_name_A", all = TRUE)
pro_score_m2<-pro_score_m
pro_score_m2[,2:4][is.na(pro_score_m2[,2:4])]<-0
pro_score_m2$total_score<-rowSums(pro_score_m2[,2:4])
pro_score_m3<-select(pro_score_m2, Gene_name_A, total_score)
colnames(pro_score_m3)<-c("Gene_name", "total_score")


# Network scoring #

pro_score_m4<-pro_score_m3
pro_score_m4$new_score<-NA


pb = txtProgressBar(min = 0, max = length(pro_score_m3$Gene_name), initial = 0)
for (i in 1:length(pro_score_m3$Gene_name)){
  scorebiogrid<-0
  scoreintact<-0
  scorewb<-0
  
  aw<-biogrid_gnames[which(biogrid_gnames$Gene_name_A %fin% pro_score_m3$Gene_name[i]),1]
  aw1<-biogrid_gnames[which(biogrid_gnames$Gene_name_B %fin% pro_score_m3$Gene_name[i]),2]
  awx<-c(aw,aw1)
  
  if (length(awx) != 0){
    scorebiogrid<-sum(na.omit(pro_score_m3$total_score[fmatch(awx, pro_score_m3$Gene_name)]))/length(awx)
  }
  
  aw2<-intact_gnames[which(intact_gnames$Gene_name_A %fin% pro_score_m3$Gene_name[i]),1]
  aw3<-intact_gnames[which(intact_gnames$Gene_name_B %fin% pro_score_m3$Gene_name[i]),2]
  awy<-c(aw2,aw3)
  
  if (length(awy) != 0){
    scoreintact<-sum(na.omit(pro_score_m3$total_score[fmatch(awy, pro_score_m3$Gene_name)]))/length(awy)
  }
  
  aw4<-wormbaseP1[which(wormbaseP1$X1 %fin% pro_score_m3$Gene_name[i]),2]
  aw5<-wormbaseP1[which(wormbaseP1$X2 %fin% pro_score_m3$Gene_name[i]),1]
  awz<-c(aw4,aw5)
  
  if (length(awz) != 0){
    scorewb<-sum(na.omit(pro_score_m3$total_score[fmatch(awz, pro_score_m3$Gene_name)]))/length(awz)
  }
  pro_score_m4$new_score[i]<-scorebiogrid + scoreintact + scorewb
  setTxtProgressBar(pb,i)
}


pro_score_m4$network_weighted<-pro_score_m4$new_score/2
pro_score_m4$Final_score<-rowSums(pro_score_m4[,c(2,4)])
colnames(wormbaseP1)[1:2]<-c("Interactor_A", "Interactor_B")


# Cluster

aa<-read.table("aa.txt", stringsAsFactors = FALSE)
aa<-gnameConverter(aa, "gene_name")
colnames(aa)<-c("Gene_name", "cluster_number")
clstr_with_scores<-merge(aa, cluster_scores, by = "cluster_number")

clstr_with_scores2<-clstr_with_scores[order(clstr_with_scores$score, decreasing = TRUE),]
clstr_with_scores2<-clstr_with_scores2[!duplicated(clstr_with_scores2$Gene_name),]
aa<-aa[-12842,]

write.table(aa, "aa.csv", row.names = FALSE, quote = FALSE, sep = "\t")


# Motifs

LM1_RFX1<-gnameConverter(LM1_RFX1, "gene_name")
RFX1<-gnameConverter(RFX1, "gene_name")
RFX1_2<-gnameConverter(RFX1_2, "gene_name")
foxj1<-gnameConverter(foxj1, "gene_name")
rfx1_transfac<-gnameConverter(rfx1_transfac, "gene_name")

motifss<-data.frame("Names" = c("LM1_RFX1", "RFX1", "RFX1_2", "rfx1_transfac"))#, "sources" = c(LM1_RFX1_pvalue, RFX1_pvalue, RFX1_2_pvalue, rfx1_transfac_pvalue))

chisqtable<-data.frame("in"=c(length(LM1_RFX1$gene_name[which(LM1_RFX1$gene_name %in% ciliaryGenes1$Gene.Name)]), length(LM1_RFX1$gene_name[which(LM1_RFX1$gene_name %in% not_ciliary$Gene.Name)])), 
                       "out" = c((436-length(LM1_RFX1$gene_name[which(LM1_RFX1$gene_name %in% ciliaryGenes1$Gene.Name)])),(972-length(LM1_RFX1$gene_name[which(LM1_RFX1$gene_name %in% not_ciliary$Gene.Name)]))))
motifss$sources[1]<-chisq.test(chisqtable, correct = F)$p.value

chisqtable<-data.frame("in"=c(length(RFX1$gene_name[which(RFX1$gene_name %in% ciliaryGenes1$Gene.Name)]), length(RFX1$gene_name[which(RFX1$gene_name %in% not_ciliary$Gene.Name)])), 
                       "out" = c((436-length(RFX1$gene_name[which(RFX1$gene_name %in% ciliaryGenes1$Gene.Name)])),(972-length(RFX1$gene_name[which(RFX1$gene_name %in% not_ciliary$Gene.Name)]))))
motifss$sources[2]<-chisq.test(chisqtable, correct = F)$p.value

chisqtable<-data.frame("in"=c(length(RFX1_2$gene_name[which(RFX1_2$gene_name %in% ciliaryGenes1$Gene.Name)]), length(RFX1_2$gene_name[which(RFX1_2$gene_name %in% not_ciliary$Gene.Name)])), 
                       "out" = c((436-length(RFX1_2$gene_name[which(RFX1_2$gene_name %in% ciliaryGenes1$Gene.Name)])),(972-length(RFX1_2$gene_name[which(RFX1_2$gene_name %in% not_ciliary$Gene.Name)]))))
motifss$sources[3]<-chisq.test(chisqtable, correct = F)$p.value

chisqtable<-data.frame("in"=c(length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% ciliaryGenes1$Gene.Name)]), length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% not_ciliary$Gene.Name)])), 
                       "out" = c((436-length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% ciliaryGenes1$Gene.Name)])),(972-length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% not_ciliary$Gene.Name)]))))
motifss$sources[4]<-chisq.test(chisqtable, correct = F)$p.value

motifss$log<-(1/log2(motifss$sources))+0.1
motifss$score<-motifss$log*10

LM1_RFX1$RFX1_score1<-motifss$score[1]
RFX1$RFX1_score2<-motifss$score[2]
RFX1_2$RFX1_score3<-motifss$score[3]
rfx1_transfac$rfx1_cur_score<-motifss$score[4]

LM1_RFX1<-unique(LM1_RFX1)
RFX1<-unique(RFX1)
RFX1_2<-unique(RFX1_2)
rfx1_transfac<-unique(rfx1_transfac)


motif_score<-merge(LM1_RFX1, RFX1, by = "gene_name", all = TRUE) %>%
  merge(RFX1_2, by = "gene_name", all = TRUE) %>%
  #merge(foxj1, by = "gene_name", all = TRUE) %>%
  merge(rfx1_transfac, by = "gene_name", all = TRUE)
motif_score<-aggregate(motif_score[,2:5], by=list(Gene_name=motif_score$gene_name), FUN= function(x) {if (length(unique(x))==1){x[1]} else {x[which(!is.na(x))]}})

motif_score[,2:5][is.na(motif_score[,2:5])]<-0
colnames(motif_score)[1]<-"Gene_name"
motif_score$total_score<-rowSums(motif_score[,2:5])


# Publication scores

all_pub2$Gene_name<-gsub("\t\t", "", all_pub2$Gene_name)
all_pub2$Publication[which(all_pub2$Pub_number == 55)]<-"Marshall et al."

all_pub2$Gene_name[which(all_pub2$Gene_name == "37257.0")]<-"SMIM4"
all_pub2$Gene_name[which(all_pub2$Gene_name == "3|KIF1B")]<-"KIF1B"
all_pub2$Pub_number[which(all_pub2$Publication == "Blacque et al._sage")]<-56


aa_pub2<-data.frame(Gene_name = unique(all_pub2$Gene_name), stringsAsFactors = FALSE)
for (i in c(1:56)){
  
  b<-all_pub2$Gene_name[which(all_pub2$Pub_number == i)]
  
  if (length(unique(all_pub2$Publication[which(all_pub2$Pub_number == i)])) != 0){
    aa_pub2[,i+1]<-0
    aa_pub2[which(aa_pub2$Gene_name %in% b),i+1]<-1
    colnames(aa_pub2)[i+1]<- unique(all_pub2$Publication[which(all_pub2$Pub_number == i)])
  }
  else {aa_pub2[,i+1]<-0}
}
aa_pub2<-aa_pub2[,c(-35,-37)]


drosophila<-fread("drosophila_homolog.txt", header = TRUE)
drosophila$`Gene Symbol`[which(drosophila$`Gene Symbol` == "")]<-NA
s<-strsplit(as.character(drosophila$`Gene Symbol`), '; ')
drosophila2<-data.frame(Gene_name = rep(drosophila$`Gene Symbol and Synonyms`, sapply(s, FUN=length)), Gene_syn=unlist(s), stringsAsFactors = FALSE)

aa_pub3<-merge(aa_pub2, drosophila2, by = "Gene_name", all = TRUE)
aa_pub3$Gene_name[which(!is.na(aa_pub3$Gene_syn))]<-aa_pub3$Gene_syn[which(!is.na(aa_pub3$Gene_syn))]
aa_pub3<-aa_pub3[,-56]
aa_pub4<-aggregate(aa_pub3[,2:55], by=list(Gene_name=aa_pub3$Gene_name), FUN=sum)

aa_pub4<-gnameConverter(aa_pub4, "Gene_name")
aa_pub4<-aggregate(aa_pub4[,2:55], by=list(Gene_name=aa_pub4$Gene_name), FUN=sum)


for (i in 2:55){
  aa_pub4[[i]][which(aa_pub4[[i]]!=0)]<-1
  aa_pub4[[i]][is.na(aa_pub4[[i]])]<-0
}

eee1<-data.frame(Paper = colnames(aa_pub4)[2:55])
eee1$score<-0
aa_pub4<-data.frame(aa_pub4, stringsAsFactors = FALSE)
for (i in 2:55){
  eee1$score[i-1]<-sum(aa_pub4[aa_pub4$Gene_name %in% ciliaryGenes1$Gene.Name, i])/length(ciliaryGenes1$Gene.Name)    #sum(aa1[,i])
}

for (i in 2:55){
  eee1$ciliary_genes[i-1]<-length(which(aa_pub4$Gene_name[aa_pub4[[i]]>0] %in% ciliaryGenes1$Gene.Name))
  eee1$not_ciliary_genes[i-1]<-length(which(aa_pub4$Gene_name[aa_pub4[[i]]>0] %in% not_ciliary$Gene.Name))
}

for (i in 1:54){
  chisqtable<-data.frame("in"=c(eee1$ciliary_genes[i], eee1$not_ciliary_genes[i]), "out" = c((436-eee1$ciliary_genes[i]),(974-eee1$not_ciliary_genes[i])))
  eee1$pvalue[i]<-chisq.test(chisqtable, correct = F)$p.value
}

for (i in 2:55){
  aa_pub4[[i]]<-aa_pub4[[i]]*eee1$score[which(eee1$Paper == colnames(aa_pub4)[i])]
}


aa_pub4$scores<-rowSums(aa_pub4[,2:55])
pub_scores<-select(aa_pub4, Gene_name, scores)


# Protein atlas score

protein_atlas_score<-fread("protein_atlas_scores.csv")
protein_atlas_score<-gnameConverter(protein_atlas_score, "Gene_name")


# single cell and genetic int. score

gen_score<-select(glist_wormbase2, Gene_name, truescore)
sc_score<-select(singleCellscore_u, Gene_name, Single_cell_score)
colnames(gen_score)<-c("Gene_name", "total_score")
colnames(sc_score)<-c("Gene_name", "total_score")



# Merge 'em all

merged_scores<-merge(pro_score_m4[,c(1,5)], gen_score, by = "Gene_name", all = TRUE)
merged_scores<-merge(merged_scores, sc_score, by = "Gene_name", all = TRUE)
colnames(merged_scores)<-c("Gene_name", "Protein_interaction_score", "Genetic_interaction_score", "Single_cell_score")
merged_scores2<-merge(merged_scores, clstr_with_scores2[,2:3], by = "Gene_name", all = TRUE)
merged_scores3<-merge(merged_scores2, motif_score[,c(1,6)], by = "Gene_name", all = TRUE)
merged_scores3<-merge(merged_scores3, pub_scores, by = "Gene_name", all = TRUE)
merged_scores3<-merge(merged_scores3, protein_atlas_score, by = "Gene_name", all = TRUE)
colnames(merged_scores3)[5:8]<-c("Cluster_score", "Motif_score", "Publication_score", "Protein_atlas_score")

merged_scores3[,c(2:8)][is.na(merged_scores3[,c(2:8)])]<-0

merged_scores3$Total_score<-rowSums(merged_scores3[,2:8])

merged_scores3$Protein_interaction_score[merged_scores3$Protein_interaction_score>=0]<-normalization(merged_scores3$Protein_interaction_score[merged_scores3$Protein_interaction_score>=0])
merged_scores3$Protein_interaction_score[merged_scores3$Protein_interaction_score<=0]<-stdize(merged_scores3$Protein_interaction_score[merged_scores3$Protein_interaction_score<=0])
merged_scores3$Single_cell_score<-normalization(merged_scores3$Single_cell_score)
merged_scores3$Cluster_score<-normalization(merged_scores3$Cluster_score)
merged_scores3$Motif_score<-normalization(merged_scores3$Motif_score)
merged_scores3$Publication_score<-normalization(merged_scores3$Publication_score)

merged_scores3$Norm_total_score<-rowSums(merged_scores3[,2:8])


# Decide weighting of each section

sections<-data.frame(section_names = colnames(merged_scores3)[2:8])
sections$score<-0

for (i in 2:4){
  sections$score[i-1]<-length(which(merged_scores3$Gene_name[merged_scores3[[i]]>0] %in% ciliaryGenes1$Gene.Name))/length(ciliaryGenes1$Gene.Name)    #sum(aa1[,i])
}
sections$score[4]<-length(which(merged_scores3$Gene_name[merged_scores3[[5]]>=0.2] %in% ciliaryGenes1$Gene.Name))/length(ciliaryGenes1$Gene.Name)    #sum(aa1[,i])
sections$score[5]<-length(which(merged_scores3$Gene_name[merged_scores3[[6]]>0] %in% ciliaryGenes1$Gene.Name))/length(ciliaryGenes1$Gene.Name)    #sum(aa1[,i])
sections$score[6]<-length(which(merged_scores3$Gene_name[merged_scores3[[7]]>0] %in% ciliaryGenes1$Gene.Name))/length(ciliaryGenes1$Gene.Name)    #sum(aa1[,i])
sections$score[7]<-length(which(merged_scores3$Gene_name[merged_scores3[[8]]>0] %in% ciliaryGenes1$Gene.Name))/length(ciliaryGenes1$Gene.Name)

norm_merged4<-merged_scores3
for (i in 2:8){
  norm_merged4[[i]]<-norm_merged4[[i]]*sections$score[i-1]
}

norm_merged4$Weighted_total_scores<-rowSums(norm_merged4[,2:8])

norm_merged5<-mousegnameConverter(norm_merged4, "Gene_name")
norm_merged5[duplicated(norm_merged5$Gene_name),]
norm_merged5<-aggregate(norm_merged5[2:11], by=list(Gene_name=norm_merged5$Gene_name), FUN = sum)
norm_merged6<-norm_merged5[order(norm_merged5$Weighted_total_scores, decreasing = TRUE),]
norm_merged6$Seq<-c(1:length(norm_merged6$Weighted_total_scores))



sections$bayes_score[1]<-length(which(merged_scores3$Gene_name[merged_scores3[[2]]>0] %in% ciliaryGenes1$Gene.Name))/length(pro_score_m4$Gene_name[which(pro_score_m4$Final_score > 0)])    #sum(aa1[,i])
sections$bayes_score[2]<-length(which(merged_scores3$Gene_name[merged_scores3[[3]]>0] %in% ciliaryGenes1$Gene.Name))/length(gen_score$Gene_name[which(gen_score$total_score > 0)])    #sum(aa1[,i])
sections$bayes_score[3]<-length(which(merged_scores3$Gene_name[merged_scores3[[4]]>0] %in% ciliaryGenes1$Gene.Name))/length(sc_score$Gene_name)    #sum(aa1[,i])
sections$bayes_score[4]<-length(which(merged_scores3$Gene_name[merged_scores3[[5]]>=0.2] %in% ciliaryGenes1$Gene.Name))/length(clstr_with_scores2$Gene_name[which(clstr_with_scores2$score>=2)])    #sum(aa1[,i])
sections$bayes_score[5]<-length(which(merged_scores3$Gene_name[merged_scores3[[6]]>0] %in% ciliaryGenes1$Gene.Name))/length(motif_score$Gene_name)    #sum(aa1[,i])
sections$bayes_score[6]<-length(which(merged_scores3$Gene_name[merged_scores3[[7]]>0] %in% ciliaryGenes1$Gene.Name))/length(pub_scores$Gene_name)    #sum(aa1[,i])
sections$bayes_score[7]<-length(which(merged_scores3$Gene_name[merged_scores3[[8]]>0] %in% ciliaryGenes1$Gene.Name))/length(protein_atlas_score$Gene_name)

norm_merged7<-norm_merged6[which(norm_merged6$Gene_name %in% gene_synonyms2$Gene_name),]
norm_merged7$Seq<-c(1:length(norm_merged7$Weighted_total_scores))

write.table(norm_merged7, "ciliogenics_ordered_list.csv", row.names = FALSE, quote = FALSE, sep = ",")


# Things for website

biogrid_gnamesW<-biogrid_gnames[,c(2,1)]
biogrid_gnamesTemp<-biogrid_gnamesW[,c(2,1)]
colnames(biogrid_gnamesTemp)<-colnames(biogrid_gnamesW)
biogrid_gnames_web<-rbind(biogrid_gnamesW, biogrid_gnamesTemp)
biogrid_gnames_web$type<-"Biogrid"
biogrid_gnames_web<-unique(biogrid_gnames_web)

intact_gnamesW<-intact_gnames[,c(2,1)]
intact_gnamesTemp<-intact_gnamesW[,c(2,1)]
colnames(intact_gnamesTemp)<-colnames(intact_gnamesW)
intact_gnames_web<-rbind(intact_gnamesW, intact_gnamesTemp)
intact_gnames_web$type<-"Intact"
intact_gnames_web<-unique(intact_gnames_web)

wormbaseP1W<-wormbaseP1[,c(2,1)]
wormbaseP1Temp<-wormbaseP1W[,c(2,1)]
colnames(wormbaseP1Temp)<-c("Gene_name_A", "Gene_name_B")
colnames(wormbaseP1W)<-c("Gene_name_A", "Gene_name_B")
wormbaseP1_web<-rbind(wormbaseP1W, wormbaseP1Temp)
wormbaseP1_web$type<-"Wormbase"
wormbaseP1_web<-unique(wormbaseP1_web)

colnames(nscores2)[1]<-"Gene_name"
nscores2<-merge(nscores2, aa, by = "Gene_name")
write.table(nscores2, "nscores21.csv", row.names = FALSE, quote = FALSE, sep = "\t")

write.table(biogrid_gnames_web, "biogrid.csv", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(intact_gnames_web, "intact.csv", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(wormbaseP1_web, "wormbaseP1.csv", row.names = FALSE, quote = FALSE, sep = "\t")



write.table(all_pub2, "publications_original.csv", row.names = FALSE, quote = FALSE, sep = "|")
all_pub3<-gnameConverter(all_pub2, "Gene_name")
write.table(all_pub3[,-3], "publications.csv", row.names = FALSE, quote = FALSE, sep = "|")

write.table(gene_synonyms2, "gene_synonyms2.csv", row.names = FALSE, quote = FALSE, sep = "\t")


# ***************************************************




for (i in c(2:4,6)){
  sections$ciliary_genes[i-1]<-length(which(merged_scores3$Gene_name[merged_scores3[[i]]>0] %in% ciliaryGenes1$Gene.Name))
  sections$not_ciliary_genes[i-1]<-length(which(merged_scores3$Gene_name[merged_scores3[[i]]>0] %in% not_ciliary$Gene.Name))
}
for (i in 5){
  sections$ciliary_genes[i-1]<-length(which(merged_scores3$Gene_name[merged_scores3[[i]]>=0.2] %in% ciliaryGenes1$Gene.Name))
  sections$not_ciliary_genes[i-1]<-length(which(merged_scores3$Gene_name[merged_scores3[[i]]>=0.2] %in% not_ciliary$Gene.Name))
}
for (i in 7){
  sections$ciliary_genes[i-1]<-length(which(merged_scores3$Gene_name[merged_scores3[[i]]>0] %in% ciliaryGenes1$Gene.Name))
  sections$not_ciliary_genes[i-1]<-length(which(merged_scores3$Gene_name[merged_scores3[[i]]>0] %in% not_ciliary$Gene.Name))
}

ciliarymerged<-merged_scores3[which(merged_scores3$Gene_name %in% ciliaryGenes1$Gene.Name),]
notciliarymerged<-merged_scores3[which(merged_scores3$Gene_name %in% not_ciliary$Gene.Name),]

for (i in 1:6){
  chisqtable<-data.frame("in"=c(sections$ciliary_genes[i], sections$not_ciliary_genes[i]), "out" = c((421-sections$ciliary_genes[i]),(933-sections$not_ciliary_genes[i])))
  sections$pvalue[i]<-chisq.test(chisqtable, correct = F)$p.value
}

sections$section_names<-as.character(sections$section_names)
write.table(sections, "stats_sections.csv", quote = FALSE, sep = ",", row.names = FALSE)

norm_merged4<-select(merged_scores3, Gene_name)

for (i in 2:length(merged_scores3)){
  norm_merged4[[i]]<-normalization(merged_scores3[[i]])
  norm_merged4[[i]]<-norm_merged4[[i]]*sections$score[i-1]
  colnames(norm_merged4)[i]<-sections$section_names[i-1]
}

norm_merged4$total_score<-rowSums(norm_merged4[,2:7])
sorted_norm_merged4<-norm_merged4[order(-total_score),]
rn_sorted_norm_merged4<-sorted_norm_merged4
rn_sorted_norm_merged4$order<-rownames(rn_sorted_norm_merged4)


