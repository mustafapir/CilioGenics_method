library(utils)
library(clipr)
library(readxl)
library(data.table)
library(stringr)
library(tidyverse)
#install.packages("pdftools")
library(pdftools)

pub<-read_xlsx("CilioGenics.xlsx")
mouse_human<-fread("mouse_human_homology_all.csv")
mouse_human<-mouse_human[,c(2,4)]
colnames(mouse_human)<-c("Gene_name", "mouse_gene_name")
celegans_human<-read_xlsx("C_elegans_Human_Homologs_20.01.2019.xlsx")


# Partial matching function

`%rin%` = function (pattern, list) {
  vapply(pattern, function (p) any(grepl(p, list)), logical(1L), USE.NAMES = FALSE)
}





wHumSymb<-pub[which(pub$`HUMAN SYMBOL` != "-"),]
wHumSymb<-wHumSymb[which(!is.na(wHumSymb$`HUMAN SYMBOL`)),]

upCase<-c()
for (i in 1:length(wHumSymb$`HUMAN SYMBOL`)){
  a<-wHumSymb$`HUMAN SYMBOL`[i]
  if (toupper(a)==a){
    upCase<-rbind(upCase,wHumSymb[i,])
  }
}

pubGenes<-data.frame(Gene_name = upCase$`HUMAN SYMBOL`, Publication = upCase$PROJECT, stringsAsFactors = FALSE)


lowCase<-c()
for (i in 1:length(wHumSymb$`HUMAN SYMBOL`)){
  a<-wHumSymb$`HUMAN SYMBOL`[i]
  if (toupper(a)!=a){
    lowCase<-rbind(lowCase,wHumSymb[i,])
  }
}

mouseSymbol<-lowCase$`HUMAN SYMBOL`
write_clip(mouseSymbol)

lowToUp<-fread("lowCaseToUpCase.txt")
lowToUp<-lowToUp[,-3]
lowToUp<-lowToUp[-1,]

#s<-strsplit(as.character(lowToUp$V2), ';')
#lowcaseUpcase<-data.frame(mouse_gene_name = rep(lowToUp$V1, sapply(s, FUN=length)), Gene_name=unlist(s))


lowCase1<-data.frame(mouse_gene_name = lowCase$`HUMAN SYMBOL`, Publication = lowCase$PROJECT, stringsAsFactors = FALSE)
lowCase2<-merge(lowCase1, mouse_human, by = "mouse_gene_name", all = TRUE)
#lowCase3<-merge(lowCase1, mouse_human, by = "mouse_gene_name")

orffromlowCase2<-lowCase2[grep("orf", lowCase2$mouse_gene_name),]
orffromlowCase2<-orffromlowCase2[c(-482,-483),]
lowCase2$Gene_name[which(lowCase2$mouse_gene_name %in% orffromlowCase2$mouse_gene_name)]<-orffromlowCase2$mouse_gene_name

lowCase3<-lowCase2[-which(is.na(lowCase2$Publication)),]
lowCase3toNotKnown<-lowCase3[which(is.na(lowCase3$Gene_name)),]
write_clip(lowCase3toNotKnown$mouse_gene_name)
lowCase3<-lowCase3[-which(is.na(lowCase3$Gene_name)),]


lowCase_Rik<-fread("lowCase-Rik.txt")
lowCase_Rik<-lowCase_Rik[,1:2]
lowCase_Rik<-lowCase_Rik[-which(lowCase_Rik$V2 == ""),]
s<-strsplit(as.character(lowCase_Rik$V2), ';')
lowCase_Rik<-data.frame(mouse_gene_name = rep(lowCase_Rik$V1, sapply(s, FUN=length)), Gene_name=unlist(s))
lowCase_Rik<-merge(lowCase1, lowCase_Rik, by = "mouse_gene_name")

lowCaseTotal<-rbind(lowCase3, lowCase_Rik)

whsymbolTotal<-rbind(pubGenes, lowCaseTotal[,2:3])
whsymbolTotal<-unique(whsymbolTotal)

write.csv(unique(whsymbolTotal$Publication), "whSymbolTotal_pub.csv", quote = FALSE, row.names = FALSE)
wHsymbol_pub_num<-fread("wHsymbol_pub_num.txt")
whsymbolTotal_wPubNo<-merge(whsymbolTotal, wHsymbol_pub_num, by = "Publication")

s<-strsplit(as.character(whsymbolTotal_wPubNo$Gene_name), '; ')
whsymbolTotal_wPubNo<-data.frame(Publication = rep(whsymbolTotal_wPubNo$Publication, sapply(s, FUN=length)), 
                                 Pub_number = rep(whsymbolTotal_wPubNo$Pub_number, sapply(s, FUN=length)), 
                                 Gene_name=unlist(s))

s<-strsplit(as.character(whsymbolTotal_wPubNo$Gene_name), ',')
whsymbolTotal_wPubNo<-data.frame(Publication = rep(whsymbolTotal_wPubNo$Publication, sapply(s, FUN=length)), 
                                 Pub_number = rep(whsymbolTotal_wPubNo$Pub_number, sapply(s, FUN=length)), 
                                 Gene_name=unlist(s))

whsymbolTotal_wPubNo$Publication<-papers$V2[match(whsymbolTotal_wPubNo$Pub_number, papers$V1)]
whsymbolTotal_final<-select(whsymbolTotal_wPubNo, Publication, Gene_name)
whsymbolTotal_final$Gene_ID<-whsymbolTotal_final$Gene_name



##********************************************************************************************************
#  Without Human Symbol #

woHumSymb<-pub[which(pub$`HUMAN SYMBOL` == "-"),]
woHumSymb1<-pub[which(is.na(pub$`HUMAN SYMBOL`)),]
woHumSymb<-rbind(woHumSymb, woHumSymb1)

NP<-woHumSymb[grep("NP_", woHumSymb$ENSEMBLE_GENE_ID),]
NP<-NP[,c(1,3)]
NP1<-separate(NP, ENSEMBLE_GENE_ID, "NP_ID", sep = " \\(E-")
ncbi_np<-fread("ncbi_np.txt") %>% select(`Protein product`, Locus)
ncbi_np<-separate(ncbi_np, `Protein product`, "NP_ID", sep = "\\.")
NP_final<-merge(NP1, ncbi_np)
names(NP_final)<-c("Gene_ID", "Publication", "Gene_name")
NP_final$Pub_number<-10
NP_final$Pub_number[which(NP_final$Publication == "Andersen 2003_Centriole proteome_Table A")]<-9



swiss<-woHumSymb[grep("SWISS-PROT:", woHumSymb$ENSEMBLE_GENE_ID),]
swiss<-swiss[,c(1,3)]
swiss$ENSEMBLE_GENE_ID<-str_replace(swiss$ENSEMBLE_GENE_ID, "SWISS-PROT:", "")
write_clip(swiss$ENSEMBLE_GENE_ID)
swissConv<-fread("swiss.tab")
swissConv<-swissConv[which(toupper(swissConv$To)==swissConv$To),]
colnames(swiss)<-c("Publication", "swissID")
colnames(swissConv)<-c("swissID", "Gene_name")
swissFinal<-merge(swiss, swissConv)
swissFinal$Gene_ID<-swissFinal$Gene_name
swissFinal$Pub_number<-9


woHumSymb<-woHumSymb[-grep("NP_", woHumSymb$ENSEMBLE_GENE_ID),]
woHumSymb<-woHumSymb[-grep("SWISS-PROT:", woHumSymb$ENSEMBLE_GENE_ID),]

woHumOrg<-data.frame(unique(woHumSymb$ORGANISM))



# Mayer et al. 2007 - 51


Mayer_2007<-woHumSymb[which(woHumSymb$PROJECT == "Mayer 2007_supptable1"), c(1, 5, 7)]
Mayer_2007_human<-Mayer_2007[which(Mayer_2007$ORGANISM == "Homo sapiens"),]
Mayer_2007_rat<-Mayer_2007[which(Mayer_2007$ORGANISM == "Rattus norvegicus"),]
Mayer_2007_mouse<-Mayer_2007[which(Mayer_2007$ORGANISM == "Mus musculus"),]
Mayer_2007_rat2<-Mayer_2007[is.na(Mayer_2007$ORGANISM),]

write_clip(Mayer_2007_human$ORIGINAL_ID)
write_clip(Mayer_2007_rat$ORIGINAL_ID)
write_clip(Mayer_2007_mouse$ORIGINAL_ID)
write_clip(Mayer_2007_rat2$ORIGINAL_ID)

rat<-fread("Mayer_2007_rat.txt", header = TRUE) %>% select(`Gene Symbol and Synonyms`, `Gene Symbol`)
rat2<-fread("Mayer_2007_rat2.txt", header = TRUE) %>% select(`Gene Symbol and Synonyms`, `Gene Symbol`)
mouse<-fread("Mayer_2007_mouse.txt", header = TRUE) %>% select(`Gene Symbol and Synonyms`, `Gene Symbol`)

Mayer_2007_animal<-rbind(rat, rat2, mouse)
colnames(Mayer_2007_animal)[1]<-"ORIGINAL_ID"
Mayer_2007<-Mayer_2007[-which(Mayer_2007$ORGANISM == "Homo sapiens"),]
Mayer_2007<-merge(Mayer_2007, Mayer_2007_animal, by = "ORIGINAL_ID")

Mayer_2007_human$`Gene Symbol`<-Mayer_2007_human$ORIGINAL_ID
Mayer_2007<-rbind(Mayer_2007, Mayer_2007_human)
Mayer_2007<-Mayer_2007[-which(Mayer_2007$`Gene Symbol` == ""),]
Mayer_2007<-data.frame(Gene_name = Mayer_2007$`Gene Symbol`, Publication = Mayer_2007$PROJECT)

s<-strsplit(as.character(Mayer_2007$Gene_name), '; ')
Mayer_2007<-data.frame(Publication = rep(Mayer_2007$Publication, sapply(s, FUN=length)), Gene_name=unlist(s)) %>% unique()

Mayer_2007$Pub_number<-51


# Andersen et al. - 9

Andersen<-woHumSymb[which(woHumSymb$PROJECT == "Andersen 2003_Centriole proteome_Table A"), c(1,3)]

Andersen_Trembl<-Andersen[grep("TREMBL:", Andersen$ENSEMBLE_GENE_ID),]
Andersen_Trembl$Gene_ID<-str_replace(Andersen_Trembl$ENSEMBLE_GENE_ID, "TREMBL:", "")


Andersen_Ensembl<-Andersen[grep("ENSEMBL:", Andersen$ENSEMBLE_GENE_ID),]
Andersen_Ensembl$Gene_ID<-str_replace(Andersen_Ensembl$ENSEMBLE_GENE_ID, "ENSEMBL:", "")


Andersen_Refseq<-Andersen[grep("REFSEQ_XP:", Andersen$ENSEMBLE_GENE_ID),]
Andersen_Refseq$Gene_ID<-str_replace(Andersen_Refseq$ENSEMBLE_GENE_ID, "REFSEQ_XP:", "")

write_clip(Andersen_Trembl$Gene_ID)
write_clip(Andersen_Ensembl$Gene_ID)
write_clip(Andersen_Refseq$Gene_ID)

Andersen_tr<-fread("Andersen_trembl.tab")
colnames(Andersen_tr)<-c("Gene_ID", "Gene_name")

Andersen_final<-merge(Andersen_Trembl, Andersen_tr, by = "Gene_ID")
colnames(Andersen_final)[2]<-"Publication"
Andersen_final<-select(Andersen_final, Publication, Gene_name)
swissFinal<-select(swissFinal, Publication, Gene_name, Gene_ID)
Andersen_final<-rbind(Andersen_final, swissFinal)
Andersen_final<-unique(Andersen_final)
Andersen_final$Pub_number<-9

# Chen et al. -12

Chen<-woHumSymb[which(woHumSymb$PROJECT == "Chen 2006"), c(1,3)] %>% separate(ENSEMBLE_GENE_ID, "Gene_ID", sep = " /REP")
write_clip(Chen$Gene_ID)
Chen$Gene_name<-celegans_human$Human_Symbol[match(Chen$Gene_ID, celegans_human$Celegans_Tag)]
Chen<-select(Chen, PROJECT, Gene_name)
colnames(Chen)[1]<-"Publication"
Chen_final<-na.omit(Chen) %>% unique()
Chen_final$Pub_number<-12


# Hodges et al. -38

Hodges<-pdf_text("Hodges_ME_et.al.PDF")

Hodges1<-Hodges[1]




# Marshall et al. 2016 -55


marshall<-read_xlsx("Marshall et al 2016, Table S4.xlsx", col_names = FALSE)

marshall1<-marshall[,1]
for (i in 2:11){
  marshall1<-rbind(marshall1, setNames(marshall[,i], names(marshall1)))
}
marshall1<-na.omit(marshall1) %>% unique()
names(marshall1)<-"mouse_gene_name"
marshall_final<-merge(marshall1, mouse_human, by = "mouse_gene_name")
marshall_final$Publication<-"Marshall et al. 2016"
names(marshall_final)[1]<-"Gene_ID"
marshall_final$Pub_number<-55

# Hodges et al. -38

hodges<-read_xlsx("hodges.xlsx")
hodges_names<-read_xlsx("hodges_gene_name.xlsx")
hnms<-colnames(hodges)
ensgalp<-c()
for (i in 1:length(hodges)){
  
  #hnms<-as.character(hodges[,i])
  e<-hodges[[i]][grep("ENSGALP", hodges[[i]])]
  ensgalp<-c(ensgalp, e)
}

#************************************************************************************************
hum_map<-read.csv("mapping_human.csv", header = FALSE)
synonyms<-hum_map[,4:5]
s<-strsplit(as.character(synonyms$V5), ',')
synonyms<-data.frame(Gene_name = rep(synonyms$V4, sapply(s, FUN=length)), Gene_synonyms=unlist(s))
synonyms$Gene_id<-hum_map$V2[match(synonyms$Gene_name, hum_map$V4)]
#************************************************************************************************

nothodges<-hodges_names$Gene_name[-which(hodges_names$Gene_name %in% NPtoGene$Locus)]
nothodges<-data.frame(nothodges, stringsAsFactors = FALSE)
write_clip(nothodges$nothodges)

hodges2<-hodges_names$Gene_name[which(hodges_names$Gene_name %in% NPtoGene$Locus)]
hodges2<-data.frame(hodges2, stringsAsFactors = FALSE)

nothodges2<-nothodges$nothodges[-which(nothodges$nothodges %in% synonyms$Gene_synonyms)]
nothodges2<-data.frame(nothodges2, stringsAsFactors = FALSE)
inhodges<-synonyms$Gene_name[match(nothodges$nothodges, synonyms$Gene_synonyms)]
inhodges<-data.frame(na.omit(inhodges))

colnames(hodges2)<-colnames(inhodges)
Hodges_hum<-rbind(hodges2, inhodges)
colnames(Hodges_hum)<-"Gene_name"
Hodges_hum$Publication<-"Hodges ME. et al."
Hodges_hum$Pub_number<-38
Hodges_hum$Gene_ID<-Hodges_hum$Gene_name



write_clip(nothodges2$nothodges2)

genesForBlast<-nothodges2

ch_ref<-fread("ncbi_chlamydomonas.csv")
nothodges2$Gene_name<-ch_ref$`Protein product`[match(nothodges2$nothodges2, ch_ref$Locus)]
nothodges21<-na.omit(nothodges2)

chlamydomonas_blast<-fread("chlamydomonas_blast.txt")
nothodges21$Hum_prot<-chlamydomonas_blast$V1[match(nothodges21$Gene_name, chlamydomonas_blast$V2)]
nothodges21<-na.omit(nothodges21)

nothodges21$Rgname<-human_np$Locus[match(nothodges21$Hum_prot, human_np$`Protein product`)]

nothodges3<-select(nothodges21, nothodges2, Rgname)
colnames(nothodges3)<-c("Gene_ID", "Gene_name")

nothodges3$Publication<-Hodges_hum$Publication[5]
nothodges3$Pub_number<-Hodges_hum$Pub_number[3]

Hodges_final <- rbind(Hodges_hum, nothodges3)

all_pub <- rbind(all_pub, Hodges_final)


# Nakachi M. et al. 2011 -39

nakachi<-read_xlsx("nakachi_supp_table9_refseq.xlsx") %>%
  unique()

nakachi_forBlast.txt<-nakachi


nakachi<-read_xlsx("nakachi_table3.xlsx")
human_np<-fread("H_sapiens_ProteinTable51.txt")
nakachi$Gene_name<-human_np$Locus[match(nakachi$Human_homolog, human_np$`Protein product`)]
nakachi<-nakachi[,c(1,3)]
nakachi$Publication<-"Nakachi M. et al."
nakachi$Pub_number<- 39
nakachi<-unique(setDT(nakachi), by = "Gene_name")



all_pub<-rbind(all_pub, nakachi)
write.table(all_pub1, file = "all_pub_10.04.20.txt", quote = FALSE)

nakachi2<-fread("nakachi_table9_biodb.txt")
nakachi2<-nakachi2[-1,]
s<-strsplit(as.character(nakachi2$V2), '; ')
nakachi2<-data.frame(Gene_ID = rep(nakachi2$V1, sapply(s, FUN=length)), Gene_name=unlist(s))
nakachi2$Publication<-"Nakachi M. et al."
nakachi2$Pub_number<- 39
all_pub<-rbind(all_pub, nakachi2)

# Lauwaet et al. 2011 -40


human_np<-fread("H_sapiens_ProteinTable51.txt")

lauwaet<-fread("Giardia_lambiae.txt")
lauwaet_filtered<-lauwaet[which(lauwaet$V4>=50),]
lauwaet_filtered_sort<-lauwaet_filtered[order(-V4),]
lauwaet_f_s_unique<-unique(setDT(lauwaet_filtered_sort), by="V2")

lauwaet_f_s_unique$V1<-human_np$Locus[match(lauwaet_f_s_unique$V1, human_np$`Protein product`)]

lauwaet_final<-select(lauwaet_f_s_unique, V1, V2)
lauwaet_final$Publication<-"Lauwaet et al. 2011"
names(lauwaet_final)[1:2]<-c("Gene_name", "Gene_ID")

lauwaet_final$Pub_number<-40

# Nogales-Cadelas et al. 2009 -46


Nogales<-fread("Nogales-Cadenas R. et al.txt", header = TRUE) %>% select(1, 2)
names(Nogales)<-c("Gene_ID", "Gene_name")
s<-strsplit(as.character(Nogales$Gene_name), '; ')
Nogales<-data.frame(Gene_ID = rep(Nogales$Gene_ID, sapply(s, FUN=length)), Gene_name = unlist(s))
Nogales$Publication<-"Nogales-Cadenas R. et al."
Nogales_final<-Nogales[-which(Nogales$Gene_name == "-"),]
Nogales_final<-unique(setDT(Nogales_final), by="Gene_name")

Nogales_final$Pub_number<-46


# Cao W. et al. 2006 -47


Cao<-fread("Cao W. et al.txt")

Cao<-Cao[-which(Cao$`Gene Symbol`==""),]
Cao<-Cao[,-3]

s<-strsplit(as.character(Cao$`Gene Symbol`), '; ')
Cao_final<-data.frame(Gene_ID = rep(Cao$`GI Number`, sapply(s, FUN=length)), Gene_name = unlist(s))
Cao_final$Publication<-"Cao W. et al."

Cao_final$Pub_number<-47

# Ostrowski et al. 2002 -20


Ostrowski_a<-woHumSymb[which(woHumSymb$ORGANISM == "?"), c(1, 3)]
write_clip(Ostrowski_T12$ENSEMBLE_GENE_ID)
Ostrowski_T12<-fread("orgnotknown_uniprot.tab")
names(Ostrowski_T12)<-c("Gene_ID", "Gene_name")

Ostrowski_T12_2<-fread("Ostrowski_wo2.txt") %>% select(V1, V2)
Ostrowski_T12_2<-Ostrowski_T12_2[-which(Ostrowski_T12_2$V2 == "-")]
Ostrowski_T12_2<-Ostrowski_T12_2[-1,]
names(Ostrowski_T12_2)<-c("Gene_ID", "Gene_name")

Ostrowski_T12_3<-fread("Ostrowski_transcID.txt") %>% select(V1, V2)
Ostrowski_T12_3<-Ostrowski_T12_3[-which(Ostrowski_T12_3$V2 == "-")]
Ostrowski_T12_3<-Ostrowski_T12_3[-1,]
names(Ostrowski_T12_3)<-c("Gene_ID", "Gene_name")

Ostrowski_1<-rbind(Ostrowski_T12, Ostrowski_T12_2, Ostrowski_T12_3)
Ostrowski_1<-unique(Ostrowski_1)
Ostrowski_1$Publication<-"Ostrowski et al."

Ostrowski_1$Pub_number<-20


Ostrowski_1_rem<-Ostrowski_a[-which(Ostrowski_a$ENSEMBLE_GENE_ID %in% Ostrowski_1$Gene_ID), ]
write_clip(Ostrowski_1_rem$ENSEMBLE_GENE_ID)



Ostrowski_S1<-pub[which(pub$PROJECT == "Ostrowski 2002_Cilium proteome_TableS1"), c(1, 3)]
write_clip(Ostrowski_S1$ENSEMBLE_GENE_ID)


# Keller et al. -33

keller<-read_xlsx("Keller_2005_Table_S2.xlsx", col_names = FALSE)
keller<-unique(keller)
NPtoGene<-fread("ncbi_protein_gene_31.12.19.txt")

keller1<-separate(keller, ...1, into = "NP", sep = "\\.")
keller1<-lapply(keller1, function(x) paste0(x, "."))
keller1<-data.frame(matrix(unlist(keller1)), stringsAsFactors = FALSE)


keller1$Gene_name<-NPtoGene$Locus[pmatch(keller1$matrix.unlist.keller1.., NPtoGene$`Protein product`)]
keller1$Publication<-"Keller et al."
keller1$Pub_number<-33
colnames(keller1)[1]<-"Gene_ID"
keller_final<-na.omit(keller1)
all_pub<-rbind(all_pub, keller_final)

# Stolc et al. -34
# blast
stolc<-collection[which(collection$V6 == 34),]
write_clip(stolc$V8)

ch_ref$`Locus tag`
ch_ref$`Locus tag`<-gsub("CHLREDRAFT_", "", ch_ref$`Locus tag`)
stolc1<-data.frame(Gene_ID = stolc$V8, stringsAsFactors = FALSE)
stolc1$Prot<-ch_ref$`Protein product`[match(stolc1$Gene_ID, ch_ref$`Locus tag`)]
stolc1<-na.omit(stolc1)
stolc1$Humprot<-chlamydomonas_blast$V1[match(stolc1$Prot, chlamydomonas_blast$V2)]
stolc1<-na.omit(stolc1)
stolc1$Gene_name<-human_np$Locus[match(stolc1$Humprot, human_np$`Protein product`)]




write.csv(allPub, file = "all_pub.csv", quote = FALSE, row.names = FALSE)


# Smith et al. -36
# ?????????


# Stubbs J.L. et al. -52

stubbs<-fread("Stubbs J.L. et al.txt", fill = TRUE, header = FALSE)
colnames(stubbs)<-"Gene_name"
stubbs$Gene_ID<-stubbs$Gene_name
stubbs$Publication<-"Stubbs J.L. et al."
stubbs$Pub_number<-52

all_pub<-rbind(all_pub, stubbs)


# Arnaiz et al. -53

#    microarray data   #

paramanot <- readAAStringSet("ptetraurelia_mac_annotation_v1.transcript.fa")
paramanot <- data.frame(paramanot, stringsAsFactors = FALSE)

fileparam<- "paramecium_fasta.txt"
cat(file = fileparam, append = FALSE, sep = "\n")
for (i in 1:length(paramanot$paramanot)){
  
  a<-paste0(">",row.names(paramanot)[i],"\n",paramanot[i,1])
  cat(a, file = fileparam, append = TRUE, sep = "\n")
}


#************


microarray_par<-fread("paramecium_microarray.txt")
microarray_par_pos<- filter(microarray_par, P.Value<=0.05)
microarray_par_pos<-filter(microarray_par_pos, t>0)

arnaiz<-select(microarray_par_pos, SEQ_ID)

arnaiz$prot<-par_ref$`Protein product`[match(arnaiz$SEQ_ID, par_ref$`Locus tag`)]
arnaiz<-na.omit(arnaiz)
arnaiz$humprot<-par$V1[match(arnaiz$prot, par$V2)]
arnaiz<-na.omit(arnaiz)
arnaiz$Gene_name<-human_np$Locus[match(arnaiz$humprot, human_np$`Protein product`)]

arnaiz<-unique(arnaiz)
write_clip(arnaiz$Gene_name)

arnaiz<-arnaiz[,c(1,4)]
colnames(arnaiz)[1]<-"Gene_ID"
arnaiz$Publication<-"Arnaiz O. et al."
arnaiz$Pub_number<-53

all_pub<-rbind(all_pub, arnaiz)



# Blacque et al. -54

ens_to_gene<-fread("ens_to_gene.txt")

blacque_xbox<-read_xls("blacque_et_al_table_S1-S4.xls", sheet = 4, skip = 3)
blacque_xbox$Gene_name<-ens_to_gene$`Gene name`[match(blacque_xbox$homolog, ens_to_gene$`Protein stable ID`)]

blacque_xbox$Gene_name2<-celegans_human$Human_Symbol[pmatch(blacque_xbox$model, celegans_human$Celegans_Tag)]

bl1<-blacque_xbox[!is.na(blacque_xbox$Gene_name2), c(1, 13)]
bl2<-blacque_xbox[!is.na(blacque_xbox$Gene_name), c(1, 12)]
colnames(bl2)[2]<-colnames(bl1)[2]

bl_xbx<-rbind(bl1, bl2)
bl_xbx<-unique(setDT(bl_xbx), by= "model")
colnames(bl_xbx)<-c("Gene_ID", "Gene_name")
bl_xbx$Publication<-"Blacque et al."
bl_xbx$Pub_number<-54


blacque_sage<-read_xls("blacque_et_al_table_S1-S4.xls", sheet = 1, skip = 3)
bl_sg<-blacque_sage[which(blacque_sage$`p-value (Cil:Pan)` < 0.05 & 
                            blacque_sage$`p-value (Cil:Mus)` < 0.05 & 
                            blacque_sage$`p-value (Cil:Gut)` < 0.05),]

bl_sg1<-bl_sg[which(bl_sg$`Cil/Pan` > 1 &
                      bl_sg$`Cil/Mus` > 1 &
                      bl_sg$`Cil/Gut` > 1),]


bl_sg1$Gene_name2<-ens_to_gene$`Gene name`[match(bl_sg1$Homolog, ens_to_gene$`Protein stable ID`)]
bl_sg1$Gene_name<-celegans_human$Human_Symbol[pmatch(bl_sg1$model, celegans_human$Celegans_Tag)]

bl1_sg<-bl_sg1[!is.na(bl_sg1$Gene_name), c(1, 20)]
bl2_sg<-bl_sg1[!is.na(bl_sg1$Gene_name2), c(1, 19)]

colnames(bl2_sg)[2]<-colnames(bl1_sg)[2]

bl_sage<-rbind(bl1_sg, bl2_sg)
bl_sage<-unique(setDT(bl_sage), by= "model")

colnames(bl_sage)<-c("Gene_ID", "Gene_name")
bl_sage$Publication<-"Blacque et al."
bl_sage$Pub_number<-54

#blacque_final<-unique(setDT(blacque_final), by = "Gene_name")

bl_sage$Publication<-"Blacque et al._sage"
bl_xbx$Publication<-"Blacque et al._xbox"

blacque_final<-rbind(bl_xbx, bl_sage)

all_pub<-rbind(all_pub, blacque_final)


#**********************

write.table(all_pub, "all_publications.txt", quote = FALSE, sep = ",", row.names = FALSE)

all_pub2<-fread("all_publications.txt", stringsAsFactors = FALSE)
all_pub2$Gene_name<-as.character(all_pub2$Gene_name)
synonyms$Gene_name<-as.character(synonyms$Gene_name)

all_pub2$Gene_name<-gsub("\t\t", "", all_pub2$Gene_name)

pb = txtProgressBar(min = 1, max = length(all_pub2$Publication), initial = 1)
for (i in 1:length(all_pub2$Publication)){
  
  if (!(all_pub2$Gene_name[i] %in% hum_map$V4) && !is.na(synonyms$Gene_name[match(all_pub2$Gene_name[i], synonyms$Gene_synonyms)])){
    
    all_pub2$Gene_name[i]<-synonyms$Gene_name[match(all_pub2$Gene_name[i], synonyms$Gene_synonyms)]
  }
  setTxtProgressBar(pb,i)
}


pb = txtProgressBar(min = 1, max = length(all_pub2$Gene_name), initial = 1)

for (i in 1:length(all_pub2$Gene_name)){
  
  if (!is.na(mouse_homology$V2[match(all_pub2$Gene_name[i], mouse_homology$V4)])){
    
    all_pub2$Gene_name[i]<-mouse_homology$V2[match(all_pub2$Gene_name[i], mouse_homology$V4)]
  }
  setTxtProgressBar(pb,i)
}

all_pub2$Gene_name[which(all_pub2$Gene_name == "37257.0")]<-"SMIM4"

all_pub2$Gene_name[which(all_pub2$Gene_name == "3|KIF1B")]<-"KIF1B"

all_pub2$Publication[which(all_pub2$Pub_number == 55)]<-"Marshall et al."


all_pub2$Pub_number[which(all_pub2$Publication == "Blacque et al._sage")]<-56
all_pub2$Gene_name<-mouse_homology$V2[match(all_pub2$Gene_name, mouse_homology$V4)]

all_pub3<-all_pub2
all_pub2$Gene_name<-mousegnameConverter(all_pub2$Gene_name)
all_pub2<-gnameConverter(all_pub2, "Gene_name")
diffr1<-setdiff(all_pub2$Gene_name, all_pub3$Gene_name)

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

aa_pub3<-merge(aa_pub2, drosophila2, by = "Gene_name", all = TRUE)

aa_pub3$Gene_name[which(!is.na(aa_pub3$Gene_syn))]<-aa_pub3$Gene_syn[which(!is.na(aa_pub3$Gene_syn))]
aa_pub3<-aa_pub3[,-57]
aa_pub4<-aggregate(aa_pub3[,2:56], by=list(Gene_name=aa_pub3$Gene_name), FUN=sum)


write.table(aa_pub4, "publication_gene_scores.csv", quote = FALSE, row.names = FALSE, sep = "\t")

aa_pub4<-fread("publication_gene_scores.csv")
aa2<-aa2[-57]

aa2<-aa_pub4

aa2$scores<-rowSums(aa2[,2:55])
pub_scores<-select(aa2, Gene_name, scores)

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
  eee1$ciliary_genes[i-1]<-length(which(aa_pub4$Gene_name[aa2[[i]]>0] %in% ciliaryGenes1$Gene.Name))
  eee1$not_ciliary_genes[i-1]<-length(which(aa_pub4$Gene_name[aa2[[i]]>0] %in% not_ciliary$Gene.Name))
}

for (i in 1:54){
  chisqtable<-data.frame("in"=c(eee1$ciliary_genes[i], eee1$not_ciliary_genes[i]), "out" = c((433-eee1$ciliary_genes[i]),(971-eee1$not_ciliary_genes[i])))
  eee1$pvalue[i]<-chisq.test(chisqtable, correct = F)$p.value
}

write.table(eee1, "publication_stats.csv", quote = FALSE, sep = ",", col.names = FALSE)

mouse_homology<-fread("mouse_human_homology_all.csv")

all_pub
colnames(aa_pub4)<-colnames(aa_pub3)

for (i in 2:55){
  
  aa_pub4[[i]]<-aa_pub4[[i]]*eee1$score[which(eee1$Paper == colnames(aa_pub4)[i])]
}



aa_pub4$scores<-rowSums(aa_pub4[,2:55])
write.table(aa_pub4, "publication_gene_Withscores.csv", quote = FALSE, row.names = FALSE, sep = "\t")

pubscores2<-select(aa_pub4, Gene_name, scores)
pubscores4<-pubscores2
pubscores4<-gnameConverter(pubscores4, "Gene_name")
pubscores3<-aggregate(pubscores4['scores'], by=pubscores4['Gene_name'], sum)

drosoph<-pubscores2$Gene_name[grepl("CG", pubscores2$Gene_name, ignore.case = T)]
drosoph<-drosoph[-c(1:5)]
drosoph<-drosoph[-c(262:290)]
write_clip(drosoph)

drosophila<-fread("drosophila_homolog.txt", header = TRUE)
drosophila$`Gene Symbol`[which(drosophila$`Gene Symbol` == "")]<-NA
s<-strsplit(as.character(drosophila$`Gene Symbol`), '; ')
drosophila2<-data.frame(Gene_name = rep(drosophila$`Gene Symbol and Synonyms`, sapply(s, FUN=length)), Gene_syn=unlist(s), stringsAsFactors = FALSE)



pubscores3<-gnameConverter(pubscores2, "Gene_name")
pubscores3<-aggregate(pubscores3['scores'], by=pubscores3['Gene_name'], sum)

write.table(pubscores3, "pubscores3_9.5.20.csv", quote = FALSE, row.names = FALSE, sep = ",")
# ALL # **************************************

Mayer_2007$Gene_ID<-Mayer_2007$Gene_name
Andersen_final$Gene_ID<-Andersen_final$Gene_name
Chen_final$Gene_ID<-Chen_final$Gene_name


allPub<-rbind(Mayer_2007, Andersen_final, Chen_final, Cao_final, lauwaet_final, marshall_final, Nogales_final, NP_final, Ostrowski_1, swissFinal, myGenes_collection_final)

allPub$Publication<-papers$V2[match(allPub$Pub_number, papers$V1)]

allPub<-allPub[order(allPub$Pub_number),]



for (i in 1:length(unique(allPub$Publication))){
  
}






# Organism: "?"

orgNotKnown<-woHumSymb[which(woHumSymb$ORGANISM == "?"), c(1, 3)]
write_clip(orgNotKnown$ENSEMBLE_GENE_ID)
orgNotKnown_uprot<-fread("orgnotknown_uniprot.tab")

Broadhead_2006_Supple_Figure1<-fread("Broadhead_2006_Supple_Figure1.txt", fill = TRUE)
Alexander_E._Ivliev_et_al<-fread("Alexander_E._Ivliev_et_al.txt")


# ******************************************

collection<-fread("collection.csv")
papers<-fread("papers.csv")
collection<-select(collection, V4, V6, V8)
collection$Publication<-papers$V2[match(collection$V6, papers$V1)]

collection_woMYGENES<-collection[-which(collection$V6 %in% whsymbolTotal_wPubNo$Pub_number),]
collection_woMYGENES<-select(collection_woMYGENES, V4, V6, Publication)
colnames(collection_woMYGENES)[1:2]<-c("Gene_name", "Pub_number")

myGenes_collection<-rbind(whsymbolTotal_wPubNo, collection_woMYGENES)
myGenes_collection<-myGenes_collection[-which(myGenes_collection$Gene_name == "--"),]
myGenes_collection$Gene_ID<-myGenes_collection$Gene_name

myGenes_collection_final<-unique(myGenes_collection)



b<-sort(a$Pub_number)






lca<-unique(lowCase2$mouse_gene_name[which(is.na(lowCase2$Gene_name))])
write_clip(lca)


pubGenes<-rbind(pubGenes, lowcaseUpcase$gene_name)






woHumSymb<-pub[which(pub$`HUMAN SYMBOL` == "-"),]
woHumSymb1<-pub[which(is.na(pub$`HUMAN SYMBOL`)),]

org_list<-unique(setDT(woHumSymb1), by= "ORGANISM")



nonEns<-pub1[which(pub1$`HUMAN ENS GID` == "-"),]
withEns<-pub1[which(pub1$`HUMAN ENS GID` != "-"),]

