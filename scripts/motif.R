
options(stringsAsFactors = FALSE)


txtFile <- paste(file.path(system.file('examples', package='RcisTarget')),
                 "hypoxiaGeneSet.txt", sep="/")
geneLists <- list(hypoxia=read.table(txtFile, stringsAsFactors=FALSE)[,1])


library(igraph)
library(networkD3)
library(stringr)
library(tidyr)

# create a dataset:
data <- tibble(
  from=c("A", "A", "B", "D", "C", "D", "E", "B", "C", "D", "K", "A", "M"),
  to=c("B", "E", "F", "A", "C", "A", "B", "Z", "A", "C", "A", "B", "K")
)

# Plot
p <- simpleNetwork(data, height="100px", width="100px")
p



LM1_RFX1<-fread("LM1_RFX1") %>% select(gene_name) %>% unique()
RFX1<-fread("RFX1_M00280") %>% select(gene_name) %>% unique()
RFX1_2<-fread("RFX1_M00281") %>% select(gene_name) %>% unique()
foxj1<-fread("foxj1.txt", sep = ",", header = FALSE)
rfx1_transfac<-fread("RFX1_TRANSFAC.txt", sep = ",", header = FALSE)
xbp1<-fread("xbp-1.tsv") %>% select("Name(s) associated with the search region that contains this motif")
colnames(xbp1)[1]<-"gene_name"
wb_sq_name<-fread("wb_sq_name.txt")
colnames(wb_sq_name)[1]<-"gene_id"

cele_homology_ref<-merge(wb_sq_name, homology_celegans, by= "gene_id")


xbp12<-strsplit(xbp1$gene_name, split = "[.]") %>% unlist() %>% data.frame(stringsAsFactors = FALSE)
xbp121<-data.frame(xbp12[seq(1, 509, 2),1], stringsAsFactors = FALSE)
xbp122<-data.frame(xbp12[seq(2, 510, 2),1], stringsAsFactors = FALSE)
xbp1221<-separate(xbp122, xbp12.seq.2..510..2...1., "g", sep = 1)


xbp12211<-cbind(xbp121, xbp1221)
xbp12211$gene_name<-apply(xbp12211[,1:2], 1, paste, collapse = ".")
xbp1<-data.frame(xbp12211$gene_name, stringsAsFactors = FALSE)
colnames(xbp1)[1]<-"gene_name"

xbp11<-data.frame(pblapply(xbp1$gene_name, function(x) cele_homology_ref[,3][match(x, cele_homology_ref$`WormBase Sequence Name`)]), stringsAsFactors = FALSE)
xbp11<-data.frame(t(xbp11))
xbp11<-na.omit(xbp11) %>% unique()

#foxj1_curated<-fread("foxj1_transfac_curated.txt", sep = ",", header = FALSE)


foxj1<-apiToGene(foxj1)
rfx1_transfac<-apiToGene(rfx1_transfac)
#foxj1_curated<-apiToGene(foxj1_curated)

LM1_RFX1<-gnameConverter(LM1_RFX1, "gene_name")
RFX1<-gnameConverter(RFX1, "gene_name")
RFX1_2<-gnameConverter(RFX1_2, "gene_name")
foxj1<-gnameConverter(foxj1, "gene_name")
rfx1_transfac<-gnameConverter(rfx1_transfac, "gene_name")



LM1_RFX1$RFX1_score1<-length(LM1_RFX1$gene_name[which(LM1_RFX1$gene_name %in% ciliaryGenes1$Gene.Name)]) / length(ciliaryGenes1$Gene.Name)
RFX1$RFX1_score2<-length(RFX1$gene_name[which(RFX1$gene_name %in% ciliaryGenes1$Gene.Name)]) / length(ciliaryGenes1$Gene.Name)
RFX1_2$RFX1_score3<-length(RFX1_2$gene_name[which(RFX1_2$gene_name %in% ciliaryGenes1$Gene.Name)]) / length(ciliaryGenes1$Gene.Name)
foxj1$foxj1_score<-length(foxj1$gene_name[which(foxj1$gene_name %in% ciliaryGenes1$Gene.Name)]) / length(ciliaryGenes1$Gene.Name)
rfx1_transfac$rfx1_cur_score<-length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% ciliaryGenes1$Gene.Name)]) / length(ciliaryGenes1$Gene.Name)

for (i in 1:6){
  chisqtable<-data.frame("in"=c(sections$ciliary_genes[i], sections$not_ciliary_genes[i]), "out" = c((436-sections$ciliary_genes[i]),(933-sections$not_ciliary_genes[i])))
  sections$pvalue[i]<-chisq.test(chisqtable, correct = F)$p.value
}
chisqtable<-data.frame("in"=c(length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% ciliaryGenes1$Gene.Name)]), length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% not_ciliary$Gene.Name)])), 
                       "out" = c((436-length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% ciliaryGenes1$Gene.Name)])),(972-length(rfx1_transfac$gene_name[which(rfx1_transfac$gene_name %in% not_ciliary$Gene.Name)]))))

rfx1_transfac_pvalue<-chisq.test(chisqtable, correct = F)$p.value

motifss<-data.frame("sources" = c(LM1_RFX1_pvalue, RFX1_pvalue, RFX1_2_pvalue, rfx1_transfac_pvalue))
motifss$normalized<-normalization(motifss$log)
motifss$scaled<-scale(motifss$log)
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






stat_motif<-motif_score[which(motif_score$Gene_name %in% ciliaryGenes1$Gene.Name),]
stat_motif2<-motif_score[which(motif_score$gene_name %in% protAtlas$gene_symbol),]
stat_motif$total<-rowSums(stat_motif[,2:6])

nostatmotif<-motif_score[which(motif_score$Gene_name %in% not_ciliary$Gene.Name),]

motif_score$total_score<-rowSums(motif_score[,2:5])
control<-data.frame("a", 0, 0, 0, 0, 0, 0)
names(control)<-names(motif_score)
motif_score<-rbind(motif_score, control)
motif_score$normalized_score<-normalization(motif_score$total_score)
known_motif<-motif_score[which(motif_score$gene_name %in% ciliaryGenes1$Gene.Name)]


ord_merged_scores3<-ord_merged_scores2
colnames(ord_merged_scores3)[1]<-"gene_name"
ord_merged_scores3<-merge(ord_merged_scores3, motif_score[,c(1,8)], by = "gene_name", all = TRUE)
ord_merged_scores3[,2:13][is.na(ord_merged_scores3[,2:13])]<-0
ord_merged_scores3$total_score3<-rowSums(ord_merged_scores3[,c(11,13)])

ord_merged_scores3<-ord_merged_scores3[order(ord_merged_scores3$total_score3, decreasing = TRUE),]
ord_merged_scores3$seq3<-seq(1, length(ord_merged_scores3$gene_name), 1)

known_not_ciliary<-ord_merged_scores3[which(ord_merged_scores3$gene_name %in% negative$`Gene Name`),]
known_ciliary<-ord_merged_scores3[which(ord_merged_scores3$gene_name %in% ciliaryGenes1$Gene.Name),]

known_ciliary_motif<-motif_score[which(motif_score$gene_name %in% ciliaryGenes1$Gene.Name),]
known_not_ciliary_motif<-motif_score[which(motif_score$gene_name %in% negative$`Gene Name`),]

