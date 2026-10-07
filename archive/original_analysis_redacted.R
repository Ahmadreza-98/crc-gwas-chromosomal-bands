### Load library

setwd('<original local working directory>')

library(quantsmooth)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(haven)
library(xlsx)
library(biomaRt)
library(sva)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
library(chromPlot)
library(clusterProfiler)
library(msigdbr)
library(KEGGREST)
library(qqman)
library(RColorBrewer)
library(rtracklayer)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)
library(DESeq2)
library(maftools)
library(pheatmap)
library(sesameData)
library(wheatmap)
library(sesame)
library(AnnotationHub)
library(car)
library(scales)
library(ggpmisc)

######## Read data

data <- read.xlsx('colorectal cancer final GWAS edited.xlsx',sheetIndex = 1)
data <- data[!duplicated(data$SNPS),]


# Major band

data$major.band<- gsub('\\..*','',data$REGION)

# Minor band
data$minor.band <- data$REGION

### Major band table (Count)
final.table.major.band <- table(data$major.band) %>%
  as.data.frame()

colnames(final.table.major.band) <- c('major.band','ni')

### Minor band table (Count)
final.table.minor.band <- table(data$minor.band) %>%
  as.data.frame()

colnames(final.table.minor.band) <- c('minor.band','ni')


########## chromosome Major band ###########


chr.band.length <- read.delim2('cytoBand.txt') %>%
  subset(.,Band!='')

chr.band.length$major.band<- paste(chr.band.length$Chr,chr.band.length$Band) %>%
  gsub('chr','',.) %>%
  gsub("\\s",'',.) %>%
  gsub('\\..*','',.)

chr.major.band.length <- subset(chr.band.length,select=c('major.band','Start','End'))

chr.major.band.length.start <- chr.major.band.length %>%
  group_by(major.band) %>%
  summarise(across(c(Start,End),min))

chr.major.band.length.end <- chr.major.band.length %>%
  group_by(major.band) %>%
  summarise(across(c(Start,End),max))

chr.major.band.length <- chr.major.band.length.end
chr.major.band.length$Start <- chr.major.band.length.start$Start
chr.major.band.length$Length <- chr.major.band.length$End-chr.major.band.length$Start


########## chromosome Minor band ###########

chr.band.length$minor.band <- paste(chr.band.length$Chr,chr.band.length$Band) %>%
  gsub('chr','',.) %>%
  gsub("\\s",'',.)

chr.minor.band.length <- subset(chr.band.length,select=c('minor.band','Start','End'))

chr.minor.band.length$Length <- chr.minor.band.length$End-chr.minor.band.length$Start

###############################
###############################
###############################

### number of polymorphisms

final.table.major.band$n <- rep(nrow(data),nrow(final.table.major.band))
final.table.minor.band$n <- rep(nrow(data),nrow(final.table.minor.band))


final.table.major.band <- merge(final.table.major.band,
                                chr.major.band.length,by='major.band')

final.table.minor.band <- merge(final.table.minor.band,
                                chr.minor.band.length,by='minor.band')

#### partial chromosome bands length

final.table.major.band$theta <- final.table.major.band$Length / 3200000000
final.table.minor.band$theta <- final.table.minor.band$Length / 3200000000


#### add df1 and df2
final.table.major.band$df1 <- 2 * (final.table.major.band$n - final.table.major.band$ni +1)
final.table.major.band$df2 <- 2 * final.table.major.band$ni

final.table.minor.band$df1 <- 2 * (final.table.minor.band$n - final.table.minor.band$ni +1)
final.table.minor.band$df2 <- 2 * final.table.minor.band$ni


##### Calculate F test

final.table.major.band$F.test <- ((1-final.table.major.band$theta)/final.table.major.band$theta) * (final.table.major.band$ni/(final.table.major.band$n - final.table.major.band$ni +1))
final.table.minor.band$F.test <- ((1-final.table.minor.band$theta)/final.table.minor.band$theta) * (final.table.minor.band$ni/(final.table.minor.band$n - final.table.minor.band$ni +1))

##### add P-value and adj-Pvalue
final.table.major.band$pvalue <- pf(final.table.major.band$F.test, final.table.major.band$df1,
                                    final.table.major.band$df2 , lower.tail = FALSE)
final.table.major.band$adj.pvalue <- final.table.major.band$pvalue * nrow(final.table.major.band)


final.table.minor.band$pvalue <- pf(final.table.minor.band$F.test, final.table.minor.band$df1,
                                    final.table.minor.band$df2 , lower.tail = FALSE)
final.table.minor.band$adj.pvalue <- final.table.minor.band$pvalue * nrow(final.table.minor.band)






##### Results


write.xlsx(subset(final.table.major.band,adj.pvalue < 0.05),'major.band.results.xlsx',row.names = F)
write.xlsx(subset(final.table.minor.band,adj.pvalue < 0.05),'minor.band.results.xlsx',row.names = F)

###############################################
###### Plot significant major bands ###########

####### Histogram

df.plot <- subset(final.table.major.band, adj.pvalue < 0.05) %>%
  subset(.,select=c('major.band', 'ni', 'Start', 'End'))


df.plot$Chrom <- gsub('[a-z].*', '', df.plot$major.band) %>%
  paste0('chr', .)

# Expanding the DataFrame
df.plot <- df.plot[rep(1:nrow(df.plot), df.plot$ni), ]

# Reset row names
rownames(df.plot) <- NULL

refGeneHg <- subset(df.plot, select=c('major.band','Chrom','Start','End'))
colnames(refGeneHg) <- c('Name','Chrom',"Start" ,"End")
refGeneHg <- subset(refGeneHg,select= c('Chrom',"Start" ,"End",'Name'))


data(hg_cytoBandIdeo)
data(hg_gap)

pdf('Histogram of GWAS polymorphism on significant major bands.pdf',width =10 ,height =25 )

chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo,
          annot1=refGeneHg,figCols=8)
dev.off()

png("Histogram of GWAS polymorphism on significant major bands.png", width = 10, height = 15, units = "in", res = 900)

chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo,
          annot1=refGeneHg,figCols=8)
dev.off()

########### Distribution

df.plot <-  subset(final.table.major.band, adj.pvalue < 0.05) %>%
  subset(., select=c('major.band','Start','End'))

df.plot$Chrom <- gsub('[a-z].*', '', df.plot$major.band) %>%
  paste0('chr', .)

df.plot$Chrom <- gsub('chr','',df.plot$Chrom)

colnames(df.plot) <- c('ID',"Start" ,"End",'Chrom')

refGeneHg <- subset(df.plot,select= c('Chrom',"Start" ,"End",'ID'))


pdf('Distribution of GWAS polymorphisms on significant major bands.pdf',width =10 ,height =25)

chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo, stat= refGeneHg, statCol="Value",
          statName="Value", noHist=TRUE, figCols=8, cex=0.7, statTyp="n",
          chrSide=c(1,1,1,1,1,1,-1,1))

dev.off()


png("Distribution of GWAS polymorphisms on significant major bands.png", width = 10, height = 15, units = "in", res = 900)

chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo, stat= refGeneHg, statCol="Value",
          statName="Value", noHist=TRUE, figCols=8, cex=0.7, statTyp="n",
          chrSide=c(1,1,1,1,1,1,-1,1))
dev.off()


###########################################################################
############### Plot total polymorphisms on major band ####################

####### Histogram
df.plot <- subset(final.table.major.band,select=c('major.band', 'ni', 'Start', 'End'))


df.plot$Chrom <- gsub('[a-z].*', '', df.plot$major.band) %>%
  paste0('chr', .)

# Expanding the DataFrame
df.plot <- df.plot[rep(1:nrow(df.plot), df.plot$ni), ]

# Reset row names
rownames(df.plot) <- NULL

refGeneHg <- subset(df.plot, select=c('major.band','Chrom','Start','End'))
colnames(refGeneHg) <- c('Name','Chrom',"Start" ,"End")
refGeneHg <- subset(refGeneHg,select= c('Chrom',"Start" ,"End",'Name'))


data(hg_cytoBandIdeo)
data(hg_gap)

pdf('Histogram of GWAS polymorphism on all major bands.pdf',width =10 ,height =25 )

chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo,
          annot1=refGeneHg,figCols=8)
dev.off()

png("Histogram of GWAS polymorphism on all major bands.png", width = 10, height = 15, units = "in", res = 900)

chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo,
          annot1=refGeneHg,figCols=8)
dev.off()

########### Distribution

df.plot <- subset(final.table.major.band, select=c('major.band','Start','End'))

df.plot$Chrom <- gsub('[a-z].*', '', df.plot$major.band) %>%
  paste0('chr', .)

df.plot$Chrom <- gsub('chr','',df.plot$Chrom)

colnames(df.plot) <- c('ID',"Start" ,"End",'Chrom')

refGeneHg <- subset(df.plot,select= c('Chrom',"Start" ,"End",'ID'))


pdf('Distribution of GWAS polymorphisms on all major bands.pdf',width =10 ,height =25)

chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo, stat= refGeneHg, statCol="Value",
          statName="Value", noHist=TRUE, figCols=8, cex=0.7, statTyp="n",
          chrSide=c(1,1,1,1,1,1,-1,1))

dev.off()


png("Distribution of GWAS polymorphisms on all major bands.png", width = 10, height = 15, units = "in", res = 900)

chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo, stat= refGeneHg, statCol="Value",
          statName="Value", noHist=TRUE, figCols=8, cex=0.7, statTyp="n",
          chrSide=c(1,1,1,1,1,1,-1,1))
dev.off()



#################
#################
#################

################# Find all GWAS variations in  significant major bands
significant_GWAS_varients_GO_major <- data.frame(matrix(NA, nrow = 0, ncol = 14))
colnames(significant_GWAS_varients_GO_major) <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", 
                              "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", 
                              "geneID", "Count", "GO", 'major.band') 

significant_GWAS_varients_KEGG_major <- data.frame(matrix(NA, nrow = 0, ncol = 15))
colnames(significant_GWAS_varients_KEGG_major) <- c("category", "subcategory", "ID", "Description", "GeneRatio", 
                                "BgRatio", "RichFactor", "FoldEnrichment", "zScore", "pvalue", 
                                "p.adjust", "qvalue", "geneID", "Count", "major.band") 

for (i in subset(final.table.major.band,adj.pvalue < 0.05)$major.band){
  
  genes_name <- c()

  for (g in c(1:nrow(data[data$major.band==i,]))){
    if(data[data$major.band==i,]$Gene.Name[g]!=0) genes_name <- append(genes_name,data[data$major.band==i,]$Gene.Name[g])
    if (data[data$major.band==i,]$Gene.Name[g]==0) genes_name <- append(genes_name,data[data$major.band==i,]$Maped.downe.gene[g])
  }
  
  genes_name <- genes_name[genes_name!=0]
  
  final_gene_name <- c()
  for (n in genes_name){
    n <- gsub(' ','',n)
    final_gene_name <- append(final_gene_name,n)
  }
  
  final_gene_name <- strsplit(final_gene_name, split = ",") %>%
    unlist() %>%
    unique()

  write.table(final_gene_name,
              paste0(i,' GWAS varients.txt'),
              row.names = F, col.names = F, sep = '\t', quote = F)
  
  ####### Enrichment
  
  gene_ids <- bitr(final_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)
  
  BP <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  MF <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  CC <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  BP_dataframe <- as.data.frame(BP)
  BP_dataframe$GO <- rep('BP', nrow(BP_dataframe)) 
  BP_dataframe$major.band <- rep(i, nrow(BP_dataframe))
  
  MF_dataframe <- as.data.frame(MF)
  MF_dataframe$GO <- rep('MF', nrow(MF_dataframe)) 
  MF_dataframe$major.band <- rep(i, nrow(MF_dataframe)) 
  
  CC_dataframe <- as.data.frame(CC)
  CC_dataframe$GO <- rep('CC', nrow(CC_dataframe)) 
  CC_dataframe$major.band <- rep(i, nrow(CC_dataframe)) 
  
  GO <- rbind(BP_dataframe,MF_dataframe)
  GO <- rbind(GO,CC_dataframe)
  write.xlsx(GO,paste0(i,' GWAS varients GO.xlsx'),row.names = F)
  
  significant_GWAS_varients_GO_major <- rbind(significant_GWAS_varients_GO_major, GO)
  
  # KEGG pathway enrichment analysis
  kegg <- enrichKEGG(
    gene          = gene_ids$ENTREZID,
    organism      = "hsa",  # Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  
  kegg_dataframe <- as.data.frame(kegg)
  kegg_dataframe$major.band <- rep(i,nrow(kegg_dataframe))
  
  write.xlsx(kegg_dataframe,paste0(i,' GWAS varients KEGG.xlsx'),row.names = F)

  significant_GWAS_varients_KEGG_major <- rbind(significant_GWAS_varients_KEGG_major, kegg_dataframe)
}



######### Find all genes in significant chromosome bands (major bands)


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


significant_major_bands <- subset(final.table.major.band,adj.pvalue < 0.05,select = c('major.band','Start','End'))
significant_major_bands$Chr <- gsub('[a-z].*','',significant_major_bands$major.band)


for (i in c(1:length(significant_major_bands$major.band))){

  genes <- getBM(
    attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
    filters = c('chromosome_name', 'start', 'end','biotype'),
    values = list(significant_major_bands$Chr[i],
                  significant_major_bands$Start[i],
                  significant_major_bands$End[i],
                  "protein_coding"),
    mart = ensembl
  )
  
  genes <- genes[!genes$hgnc_symbol=='',]
  symbol <- gsub('-.*','',genes$hgnc_symbol) %>%
    unique()
  
  write.table(symbol,paste0(significant_major_bands$major.band[i],' genes.txt'),sep = '\t',quote = F,
              row.names = F,col.names = F)

}


####### Enrichment by clusterProfiler

significant_GO_major <- data.frame(matrix(NA, nrow = 0, ncol = 14))
colnames(significant_GO_major) <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", 
                              "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", 
                              "geneID", "Count", "GO", 'major.band') 

significant_KEGG_major <- data.frame(matrix(NA, nrow = 0, ncol = 15))
colnames(significant_KEGG_major) <- c("category", "subcategory", "ID", "Description", "GeneRatio", 
                               "BgRatio", "RichFactor", "FoldEnrichment", "zScore", "pvalue", 
                               "p.adjust", "qvalue", "geneID", "Count", "major.band") 


for (i in c(1:length(significant_major_bands$major.band))){
  
  genes <- getBM(
    attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
    filters = c('chromosome_name', 'start', 'end','biotype'),
    values = list(significant_major_bands$Chr[i],
                  significant_major_bands$Start[i],
                  significant_major_bands$End[i],
                  "protein_coding"),
    mart = ensembl
  )
  
  genes <- genes[!genes$hgnc_symbol=='',]
  symbol <- gsub('-.*','',genes$hgnc_symbol) %>%
    unique()
    
  gene_ids <- bitr(symbol, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)
  
  BP <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  MF <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  CC <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  BP_dataframe <- as.data.frame(BP)
  BP_dataframe$GO <- rep('BP', nrow(BP_dataframe)) 
  BP_dataframe$major.band <- rep(significant_major_bands$major.band[i], nrow(BP_dataframe))
  
  MF_dataframe <- as.data.frame(MF)
  MF_dataframe$GO <- rep('MF', nrow(MF_dataframe)) 
  MF_dataframe$major.band <- rep(significant_major_bands$major.band[i], nrow(MF_dataframe)) 
  
  CC_dataframe <- as.data.frame(CC)
  CC_dataframe$GO <- rep('CC', nrow(CC_dataframe)) 
  CC_dataframe$major.band <- rep(significant_major_bands$major.band[i], nrow(CC_dataframe)) 
  
  GO <- rbind(BP_dataframe,MF_dataframe)
  GO <- rbind(GO,CC_dataframe)
  write.xlsx(GO,paste0(significant_major_bands$major.band[i],' GO.xlsx'),row.names = F)
  significant_GO_major <- rbind(significant_GO_major,GO)
  
  # KEGG pathway enrichment analysis
  kegg <- enrichKEGG(
    gene          = gene_ids$ENTREZID,
    organism      = "hsa",  # Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  
  kegg_dataframe <- as.data.frame(kegg)
  kegg_dataframe$major.band <- rep(significant_major_bands$major.band[i],nrow(kegg_dataframe))
  
  write.xlsx(kegg_dataframe,paste0(significant_major_bands$major.band[i],' KEGG.xlsx'),row.names = F)
  significant_KEGG_major <- rbind(significant_KEGG_major,kegg_dataframe)
}

################### 
###################
###################
################### Enrich Chr major bands for significant kegg pathway genes

mutual_significant_major_band <- data.frame(matrix(NA, nrow = 0, ncol = 16))
colnames(mutual_significant_major_band) <- c("ID", "Description","size.overlap.term", "size.overlap.category", 
                                             "size.term", "size.category", "RichFactor", "FoldEnrichment", 
                                             "zScore", "pvalue", "p.adjust", "qvalue", 
                                             "geneID", "Count", "k.k","KEGG.ID")

for (i in c(1:nrow(significant_KEGG_major))){
  
  genes_in_kegg <- keggLink("hsa",significant_KEGG_major$ID[i])
  genes_in_kegg_ids <- gsub("hsa:", "", genes_in_kegg) %>% 
    as.vector()
  
  
  cp.data <- msigdbr(species= 'Homo sapiens', category='C1' )
  cp.data.subset <- select(cp.data, gs_name, entrez_gene)
  
  
  enrich.cp <- enricher(gene=genes_in_kegg_ids ,TERM2GENE = cp.data.subset )
  
  
  enrich.cp.df <- enrich.cp@result %>%
    separate(BgRatio, into= c('size.term','size.category'), sep='/') %>%
    separate(GeneRatio, into = c('size.overlap.term', 'size.overlap.category') ,sep = '/') %>%
    mutate_at(vars('size.term','size.category','size.overlap.term', 'size.overlap.category'),as.numeric) %>%
    mutate('k.k'=size.overlap.term/size.term) %>%
    .[order(.$p.adjust, decreasing = FALSE),]
  
  enrich.cp.df$Description <- gsub('chr','',enrich.cp.df$Description)
  enrich.cp.df <- enrich.cp.df[enrich.cp.df$p.adjust < 0.05,]
  
  if (as.vector(significant_KEGG_major$major.band[i]) %in% enrich.cp.df$Description){
    enrich.cp.df <- enrich.cp.df[enrich.cp.df$Description== as.vector(significant_KEGG_major$major.band[i]),]
    enrich.cp.df$"KEGG.ID" <- significant_KEGG_major$ID[i]
    mutual_significant_major_band <-rbind(mutual_significant_major_band, enrich.cp.df)
  }
  
}

mutual_significant_major_band_KEGG <- subset(mutual_significant_major_band,
                                             select= c("KEGG.ID", "Description","size.overlap.term", "size.overlap.category", 
                                                       "size.term", "size.category", "RichFactor", "FoldEnrichment", 
                                                       "zScore", "pvalue", "p.adjust", "qvalue", 
                                                       "geneID", "Count", "k.k") )

row.names(mutual_significant_major_band_KEGG) <- 1:nrow(mutual_significant_major_band_KEGG)

##### 
###################
###################
################### Enrich Chr major bands for significant GO genes

mutual_significant_major_band <- data.frame(matrix(NA, nrow = 0, ncol = 17))
colnames(mutual_significant_major_band) <- c("ID", "Description","size.overlap.term", "size.overlap.category", 
                                             "size.term", "size.category", "RichFactor", "FoldEnrichment", 
                                             "zScore", "pvalue", "p.adjust", "qvalue", 
                                             "geneID", "Count", "k.k","GO.ID","GO.Type")

for (i in c(1:nrow(significant_GO_major))){
  
  genes_in_GO <- AnnotationDbi::select(org.Hs.eg.db, keys = significant_GO_major$ID,
                                       columns = "ENTREZID", keytype = "GOALL")
  genes_in_GO <- genes_in_GO[!duplicated(genes_in_GO$ENTREZID), ]
  
  cp.data <- msigdbr(species= 'Homo sapiens', category='C1' )
  cp.data.subset <- select(cp.data, gs_name, entrez_gene)
  
  
  enrich.cp <- enricher(gene = genes_in_GO$ENTREZID ,TERM2GENE = cp.data.subset )
  
  
  enrich.cp.df <- enrich.cp@result %>%
    separate(BgRatio, into= c('size.term','size.category'), sep='/') %>%
    separate(GeneRatio, into = c('size.overlap.term', 'size.overlap.category') ,sep = '/') %>%
    mutate_at(vars('size.term','size.category','size.overlap.term', 'size.overlap.category'),as.numeric) %>%
    mutate('k.k'=size.overlap.term/size.term) %>%
    .[order(.$p.adjust, decreasing = FALSE),]
  
  enrich.cp.df$Description <- gsub('chr','',enrich.cp.df$Description)
  enrich.cp.df <- enrich.cp.df[enrich.cp.df$p.adjust < 0.05,]
  
  if (as.vector(significant_GO_major$major.band[i]) %in% enrich.cp.df$Description){
    enrich.cp.df <- enrich.cp.df[enrich.cp.df$Description== as.vector(significant_GO_major$major.band[i]),]
    enrich.cp.df$"GO.ID" <- significant_GO_major$ID[i]
    enrich.cp.df$"GO.Type" <- significant_GO_major$GO[i]
    mutual_significant_major_band <-rbind(mutual_significant_major_band, enrich.cp.df)
  }
  
}

mutual_significant_major_band_GO <- subset(mutual_significant_major_band,
                                           select= c("GO.ID","GO.Type" ,"Description","size.overlap.term", "size.overlap.category", 
                                                     "size.term", "size.category", "RichFactor", "FoldEnrichment", 
                                                     "zScore", "pvalue", "p.adjust", "qvalue", 
                                                     "geneID", "Count", "k.k") )

row.names(mutual_significant_major_band_GO) <- 1:nrow(mutual_significant_major_band_GO)
mutual_significant_major_band_GO <- mutual_significant_major_band_GO[!duplicated(mutual_significant_major_band_GO),]

#####################################################################
################### Intersect GO and KEGG ###########################

union(mutual_significant_major_band_GO$Description,
          mutual_significant_major_band_KEGG$Description)

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
################### Find all genes in significant chromosome bands (minor bands)#######################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Modify to use minor.band instead of major.band
significant_minor_bands <- subset(final.table.minor.band, adj.pvalue < 0.05, select = c('minor.band', 'Start', 'End'))
significant_minor_bands$Chr <- gsub('[a-z].*','',significant_minor_bands$minor.band)

for (i in 1:length(significant_minor_bands$minor.band)) {
  
  genes <- getBM(
    attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
    filters = c('chromosome_name', 'start', 'end', 'biotype'),
    values = list(significant_minor_bands$Chr[i],
                  significant_minor_bands$Start[i],
                  significant_minor_bands$End[i],
                  "protein_coding"),
    mart = ensembl
  )
  
  genes <- genes[!genes$hgnc_symbol == '', ]
  symbol <- gsub('-.*','', genes$hgnc_symbol) %>%
    unique()
  
  write.table(symbol, paste0(significant_minor_bands$minor.band[i], ' genes.txt'), sep = '\t', quote = F,
              row.names = F, col.names = F)
}


####### Enrichment by clusterProfiler

significant_GO_minor <- data.frame(matrix(NA, nrow = 0, ncol = 14))
colnames(significant_GO_minor) <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", 
                              "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", 
                              "geneID", "Count", "GO", 'minor.band') 

significant_KEGG_minor <- data.frame(matrix(NA, nrow = 0, ncol = 15))
colnames(significant_KEGG_minor) <- c("category", "subcategory", "ID", "Description", "GeneRatio", 
                                "BgRatio", "RichFactor", "FoldEnrichment", "zScore", "pvalue", 
                                "p.adjust", "qvalue", "geneID", "Count", "minor.band") 

for (i in 1:length(significant_minor_bands$minor.band)) {
  
  genes <- getBM(
    attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
    filters = c('chromosome_name', 'start', 'end', 'biotype'),
    values = list(significant_minor_bands$Chr[i],
                  significant_minor_bands$Start[i],
                  significant_minor_bands$End[i],
                  "protein_coding"),
    mart = ensembl
  )
  
  genes <- genes[!genes$hgnc_symbol == '', ]
  symbol <- gsub('-.*', '', genes$hgnc_symbol) %>%
    unique()
  
  gene_ids <- bitr(symbol, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)
  
  BP <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  MF <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  CC <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  BP_dataframe <- as.data.frame(BP)
  BP_dataframe$GO <- rep('BP', nrow(BP_dataframe)) 
  BP_dataframe$minor.band <- rep(significant_minor_bands$minor.band[i], nrow(BP_dataframe))
  
  MF_dataframe <- as.data.frame(MF)
  MF_dataframe$GO <- rep('MF', nrow(MF_dataframe)) 
  MF_dataframe$minor.band <- rep(significant_minor_bands$minor.band[i], nrow(MF_dataframe)) 
  
  CC_dataframe <- as.data.frame(CC)
  CC_dataframe$GO <- rep('CC', nrow(CC_dataframe)) 
  CC_dataframe$minor.band <- rep(significant_minor_bands$minor.band[i], nrow(CC_dataframe)) 
  
  GO <- rbind(BP_dataframe, MF_dataframe)
  GO <- rbind(GO, CC_dataframe)
  write.xlsx(GO, paste0(significant_minor_bands$minor.band[i], ' GO.xlsx'), row.names = F)
  significant_GO_minor <- rbind(significant_GO_minor, GO)
  
  # KEGG pathway enrichment analysis
  kegg <- enrichKEGG(
    gene          = gene_ids$ENTREZID,
    organism      = "hsa",  # Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  
  kegg_dataframe <- as.data.frame(kegg)
  kegg_dataframe$minor.band <- rep(significant_minor_bands$minor.band[i], nrow(kegg_dataframe))
  
  write.xlsx(kegg_dataframe, paste0(significant_minor_bands$minor.band[i], ' KEGG.xlsx'), row.names = F)
  significant_KEGG_minor <- rbind(significant_KEGG_minor, kegg_dataframe)
}




##################################################################
##################################################################
#################### Obtain Somatic Mutation ##################### 
##################################################################
##################################################################

## get list of project

GDC_project <- getGDCprojects()



### Mutation data
TCGA_COAD_mutation <-  GDCquery( project = "TCGA-COAD",
                                 data.category = "Simple Nucleotide Variation",
                                 data.type = "Masked Somatic Mutation",
                                 access = 'open')


TCGA_READ_mutation <- GDCquery( project = "TCGA-READ",
                                data.category = "Simple Nucleotide Variation",
                                data.type = "Masked Somatic Mutation",
                                access = 'open')
## get output

output_TCGA_COAD_mutation <- getResults(TCGA_COAD_mutation)
output_TCGA_READ_mutation <- getResults(TCGA_READ_mutation)

## download data

GDCdownload(TCGA_COAD_mutation,method = "api", directory = "maf")
GDCdownload(TCGA_READ_mutation,method = "api", directory = "maf")


TCGA_COAD_mutation_data <- GDCprepare(TCGA_COAD_mutation, directory = "maf",
                                      summarizedExperiment = F)
TCGA_READ_mutation_data <- GDCprepare(TCGA_READ_mutation, directory = "maf",
                                      summarizedExperiment = F)


##merge COAD and READ

#levels(TCGA_READ_mutation_data$X1) <- as.character(c(1:161)+457)


TCGA_CRAD_mutation_data <- rbind(TCGA_COAD_mutation_data,
                                 TCGA_READ_mutation_data)


#### Prepare maf file

TCGA_CRAD_mutation_maf <- TCGA_CRAD_mutation_data %>%
  maftools::read.maf()


####### Gene summary ########

################# get transcript length #########################

mart <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl") 

# First, map HGNC symbols to Ensembl IDs
gene_ids <- getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id"),
  filters = "hgnc_symbol",
  values = unique(TCGA_CRAD_mutation_data$Hugo_Symbol),
  mart = mart
)

# Now, fetch cds_length using ensembl_transcript_id
cds_lengths <- getBM(
  attributes = c("ensembl_transcript_id", "cds_length"),
  filters = "ensembl_transcript_id",
  values = gene_ids$ensembl_transcript_id,
  mart = mart
)

# Merge results to link HGNC symbols back
transcript_lengths <- merge(gene_ids, cds_lengths, by = "ensembl_transcript_id")
transcript_lengths <- subset(transcript_lengths,select = c('hgnc_symbol','cds_length')) %>%
  na.omit()

longest_transcripts <- transcript_lengths %>%
  group_by(hgnc_symbol) %>%
  slice_max(order_by = cds_length, n = 1) 


####  get gene summary
crad_gene_mutation <- getGeneSummary(TCGA_CRAD_mutation_maf)
colnames(longest_transcripts) <- c('Hugo_Symbol','cds_length')


### adjust

crad_gene_mutation <- merge(crad_gene_mutation,
                            longest_transcripts,
                            by= 'Hugo_Symbol')


quantiles <- quantile(crad_gene_mutation$cds_length,
                      probs = c(0.25, 0.5, 0.75))


crad_gene_mutation$total.adj <- (crad_gene_mutation$total / (crad_gene_mutation$cds_length + quantiles[[1]])) * log2(crad_gene_mutation$cds_length)
crad_gene_mutation <- crad_gene_mutation[!duplicated(crad_gene_mutation),]


###############################################
###############################################
####### Find 250 up and 250 down genes ########
###############################################
###############################################

COAD_up <- read.delim2('COAD-Adeno-Up-reg-top250.txt')
COAD_down <- read.delim2('COAD-Adeno-Down-reg-top250.txt')
COAD_all <- rbind(COAD_up,
                  COAD_down)


#### Check DEGs and mutation for significant bands


signifcant_bands_mutations_degs <- data.frame(mean_mutation = NA,
                                              count_deg = NA,
                                              band_genes= NA,
                                              band_name= NA)

for (i in c(1:length(significant_major_bands$major.band))){
  
  significant__major_bands_genes <- read.table(paste0(significant_major_bands$major.band[i],' genes.txt'),
                                               header = F)
  
  signifcant_bands_mutations_degs[nrow(signifcant_bands_mutations_degs) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% significant__major_bands_genes$V1,]$total.adj),
                                                                                    length(intersect(significant__major_bands_genes$V1,COAD_all$Gene)),
                                                                                    length(significant__major_bands_genes$V1),
                                                                                    as.character(significant_major_bands$major.band[i]))
  }


signifcant_bands_mutations_degs  <- na.omit(signifcant_bands_mutations_degs )
signifcant_bands_mutations_degs$count_deg <- as.numeric(signifcant_bands_mutations_degs$count_deg)
signifcant_bands_mutations_degs$DEGs_percent <- round(signifcant_bands_mutations_degs$count_deg/as.numeric(signifcant_bands_mutations_degs$band_genes),4)*100
signifcant_bands_mutations_degs$mean_mutation <- as.numeric(signifcant_bands_mutations_degs$mean_mutation)%>%
  round(.,5)



signifcant_bands_mutations_degs <- subset(signifcant_bands_mutations_degs, select= c("band_name","band_genes","count_deg",
                                                                                     "DEGs_percent","mean_mutation"))


signifcant_bands_mutation_data <- list()

for (i in c(1:length(significant_major_bands$major.band))){
  
  significant__major_bands_genes <- read.table(paste0(significant_major_bands$major.band[i],' genes.txt'),
                                               header = F)
  band <- as.character(significant_major_bands$major.band[i])
  signifcant_bands_mutation_data [[band]]<- crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% significant__major_bands_genes$V1,]$total.adj
}

# Create a dataframe from all groups
signifcant_bands_mutation_data <- data.frame(
  value = c(signifcant_bands_mutation_data$`10p14`,
            signifcant_bands_mutation_data$`10q25`,
            signifcant_bands_mutation_data$`11q12`,
            signifcant_bands_mutation_data$`12p13`,
            signifcant_bands_mutation_data$`15q13`,
            signifcant_bands_mutation_data$`18q21`,
            signifcant_bands_mutation_data$`19q13`,
            signifcant_bands_mutation_data$`1q41`,
            signifcant_bands_mutation_data$`20p12`,
            signifcant_bands_mutation_data$`20q13`,
            signifcant_bands_mutation_data$`6p21`,
            signifcant_bands_mutation_data$`8q24`,
            signifcant_bands_mutation_data$`9q34`),
  
  group = rep(c("10p14", 
                "10q25",
                "11q12",
                "12p13",
                "15q13",
                "18q21",
                "19q13",
                "1q41",
                "20p12",
                "20q13",
                "6p21",
                "8q24",
                "9q34"), 
              
              c(length(signifcant_bands_mutation_data$`10p14`),
                length(signifcant_bands_mutation_data$`10q25`),
                length(signifcant_bands_mutation_data$`11q12`),
                length(signifcant_bands_mutation_data$`12p13`),
                length(signifcant_bands_mutation_data$`15q13`),
                length(signifcant_bands_mutation_data$`18q21`),
                length(signifcant_bands_mutation_data$`19q13`),
                length(signifcant_bands_mutation_data$`1q41`),
                length(signifcant_bands_mutation_data$`20p12`),
                length(signifcant_bands_mutation_data$`20q13`),
                length(signifcant_bands_mutation_data$`6p21`),
                length(signifcant_bands_mutation_data$`8q24`),
                length(signifcant_bands_mutation_data$`9q34`)))
)

# T-test for finding mutation difference

mutation_t_test_results <- data.frame(band= NA,
                                      f= NA,
                                      f_test_p_value= NA,
                                      t= NA,
                                      df= NA,
                                      band_mean= NA,
                                      background_mean= NA,
                                      P_value= NA)

for (i in c(1:length(significant_major_bands$major.band))){
  
  significant__major_bands_genes <- read.table(paste0(significant_major_bands$major.band[i],' genes.txt'),
                                               header = F)
  
  background_data <- crad_gene_mutation[!crad_gene_mutation$Hugo_Symbol %in% significant__major_bands_genes$V1,]
  band_data <- crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% significant__major_bands_genes$V1,]
    
  if(var.test(band_data$total.adj, background_data$total.adj)[[3]] > 0.05) {
    
    # Variances are equal → Use var.equal = TRUE
  
    t_test_result <- t.test(band_data$total.adj, background_data$total.adj,
                            alternative = "greater", var.equal = T)
  } else {
    # Variances are unequal → Use var.equal = FALSE
    t_test_result <- t.test(band_data$total.adj, background_data$total.adj,
                            alternative = "greater", var.equal = F)
  }
  
  mutation_t_test_results[nrow(mutation_t_test_results)+1,] <- c(as.character(significant_major_bands$major.band[i]),
                                                                 var.test(band_data$total.adj, background_data$total.adj)[[1]][[1]],
                                                                 var.test(band_data$total.adj, background_data$total.adj)[[3]][[1]],
                                                                 t_test_result$statistic[[1]],
                                                                 t_test_result$parameter[[1]],
                                                                 t_test_result$estimate[[1]],
                                                                 t_test_result$estimate[[2]],
                                                                 t_test_result$p.value[[1]])
}

mutation_t_test_results <- na.omit(mutation_t_test_results)
mutation_t_test_results$P_value <- as.numeric(mutation_t_test_results$P_value)
mutation_t_test_results$adj_p_value <- p.adjust(mutation_t_test_results$P_value, method = "BH")


##### Hypergeometric test for mutation

mutation_hyper_results <- data.frame(band= NA,
                                     bands_gene= NA,
                                     bands_gene_in_mutataion_data= NA,
                                     high_mutation= NA,
                                     hyper_p_value= NA)

high_crad_gene_mutation <- crad_gene_mutation[crad_gene_mutation$total.adj >= quantile(crad_gene_mutation$total.adj)[[3]],]

for (i in c(1:length(significant_major_bands$major.band))){
  
  significant_major_bands_genes <- read.table(paste0(significant_major_bands$major.band[i],' genes.txt'),
                                               header = F)
  
  hyper_p_value <- phyper(q = length(intersect(significant_major_bands_genes$V1, high_crad_gene_mutation$Hugo_Symbol)) - 1,      
                          m = length(intersect(significant_major_bands_genes$V1,crad_gene_mutation$Hugo_Symbol)),           
                          n = length(crad_gene_mutation$Hugo_Symbol) - length(intersect(significant_major_bands_genes$V1,crad_gene_mutation$Hugo_Symbol)), 
                          k = length(high_crad_gene_mutation$Hugo_Symbol),          
                          lower.tail = FALSE) 
  
  mutation_hyper_results[nrow(mutation_hyper_results)+1,] <- c(as.character(significant_major_bands$major.band[i]),
                                                       length(significant_major_bands_genes$V1),
                                                       length(intersect(significant_major_bands_genes$V1,crad_gene_mutation$Hugo_Symbol)),
                                                       length(intersect(significant_major_bands_genes$V1, high_crad_gene_mutation$Hugo_Symbol)),
                                                       hyper_p_value)
}

mutation_hyper_results <- na.omit(mutation_hyper_results)
mutation_hyper_results$hyper_p_value <- as.numeric(mutation_hyper_results$hyper_p_value)
mutation_hyper_results$adj_p_value <- p.adjust(mutation_hyper_results$hyper_p_value, method = "BH")


# Define a vector of 12 colors
colors <- c("skyblue", "lightgreen", "yellow", "pink", "lightblue", "orange", "purple", 
            "red", "blue", "green", "brown", "gray", "cyan")

legend_names <- c("10p14", 
                  "10q25",
                  "11q12",
                  "12p13",
                  "15q13",
                  "18q21",
                  "19q13",
                  "1q41",
                  "20p12",
                  "20q13",
                  "6p21",
                  "8q24",
                  "9q34")

# Create the box plot
pdf('box plot of bands mutation.pdf')


ggplot(signifcant_bands_mutation_data, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = colors, labels = legend_names) + 
  theme_minimal() +
  labs(title = "Distribution of Adjusted Mutation Counts in Significant bands", y = "Value", fill = "Legend") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())


dev.off()



#### PNG
# Your plot code

plot <- ggplot(signifcant_bands_mutation_data, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = colors, labels = legend_names) + 
  theme_minimal() +
  labs(title = "Distribution of Adjusted Mutation Counts in Significant bands", y = "Value", fill = "Legend") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())


# Save the plot as a PNG file
ggsave("box plot of bands mutation.png", plot = plot, width = 8, height = 6, dpi = 300)



#### Hypergeometric test
degs_hyper_results <- data.frame(band= NA,
                                 bands_gene= NA,
                                 DEGs= NA,
                                 hyper_p_value= NA)
  
  
for (i in c(1:length(significant_major_bands$major.band))){
  
  significant__major_bands_genes <- read.table(paste0(significant_major_bands$major.band[i],' genes.txt'),
                                               header = F)
  
  hyper_p_value <- phyper(q = length(intersect(significant__major_bands_genes$V1,COAD_all$Gene)) - 1,      
       m = length(significant__major_bands_genes$V1),           
       n = 25000 - length(significant__major_bands_genes$V1), 
       k = length(COAD_all$Gene),          
       lower.tail = FALSE) 

  degs_hyper_results[nrow(degs_hyper_results)+1,] <- c(as.character(significant_major_bands$major.band[i]),
                                                       length(significant__major_bands_genes$V1),
                                                       length(intersect(significant__major_bands_genes$V1,COAD_all$Gene)),
                                                       hyper_p_value)
}


degs_hyper_results <- na.omit(degs_hyper_results)
degs_hyper_results$hyper_p_value <- as.numeric(degs_hyper_results$hyper_p_value)
degs_hyper_results$adj_p_value <- p.adjust(degs_hyper_results$hyper_p_value, method = "BH")

########################################################################
########################################################################



for(i in seq_along(significant_major_bands)){
  
  data[data$major.band == significant_major_bands$major.band[i], ]$Type
  
  
}


##########
##########
##########

# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", "Missense", 'UTRs (5 prime UTR and 3 prime UTR)', "Others (Synonymous, Splice Acceptor, and Stop Gained)",
               "Non-Coding", "Non-coding (annotated)", "Non-coding (uncharacterized)"),
  Count = c(1346, 632, 577, 20, 24, 11, 714, 211, 503),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub", "Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[7,'Percentage'] <- round(plot_data[7,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic", "Missense", 
                                        "Others (Synonymous, Splice Acceptor, and Stop Gained)", 
                                        "UTRs (5 prime UTR and 3 prime UTR)", "Non-coding (annotated)", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_three_levels_updated_color_v5.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
  #+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_three_levels_updated_color_v5_high_quality.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
  #+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()


##############################
##############################
##############################
##############################

# Example dot plot for KEGG

top_significant_KEGG_major <- significant_KEGG_major[order(significant_KEGG_major$p.adjust,decreasing = F),] %>%
  .[1:10,]


pdf('Top 10 KEGG Pathway Enrichment Dot Plot.pdf')

ggplot(top_significant_KEGG_major, aes(x = reorder(Description, -p.adjust), y = p.adjust, size = Count, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 KEGG Pathway Enrichment Dot Plot", x = "Pathway", y = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

dev.off()


png('Top 10 KEGG Pathway Enrichment Dot Plot.png', width = 7, height = 5, units = "in", res = 300)

ggplot(top_significant_KEGG_major, aes(x = reorder(Description, -p.adjust), y = p.adjust, size = Count, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 KEGG Pathway Enrichment Dot Plot", x = "Pathway", y = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

dev.off()



# Example dot plot for GO
# BP Plot
top_significant_BP_major <- significant_GO_major[significant_GO_major$GO =='BP',] %>%
  .[order(.$p.adjust),] %>%
  .[1:10,]

# Save as PDF
pdf('Top 10 BP Enrichment Dot Plot.pdf', width = 9, height = 5)

ggplot(top_significant_BP_major, aes(x = reorder(Description, -p.adjust), y = p.adjust, size = Count, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 Biological Process (BP) Enrichment Dot Plot", x = "Biological Process", y = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

dev.off()

# Save as PNG
png('Top 10 BP Enrichment Dot Plot.png', width = 9, height = 5, units = "in", res = 300)

ggplot(top_significant_BP_major, aes(x = reorder(Description, -p.adjust), y = p.adjust, size = Count, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 Biological Process (BP) Enrichment Dot Plot", x = "Biological Process", y = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

dev.off()

# MF Plot
top_significant_MF_major <- significant_GO_major[significant_GO_major$GO =='MF',] %>%
  .[order(.$p.adjust),] %>%
  .[1:10,]

# Save as PDF
pdf('Top 10 MF Enrichment Dot Plot.pdf', width = 8, height = 5)

ggplot(top_significant_MF_major, aes(x = reorder(Description, -p.adjust), y = p.adjust, size = Count, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 Molecular Function (MF) Enrichment Dot Plot", x = "Molecular Function", y = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

dev.off()

# Save as PNG
png('Top 10 MF Enrichment Dot Plot.png', width = 8, height = 5, units = "in", res = 300)

ggplot(top_significant_MF_major, aes(x = reorder(Description, -p.adjust), y = p.adjust, size = Count, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 Molecular Function (MF) Enrichment Dot Plot", x = "Molecular Function", y = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

dev.off()

# CC Plot
top_significant_CC_major <- significant_GO_major[significant_GO_major$GO =='CC',] %>%
  .[order(.$p.adjust),] %>%
  .[1:10,]

# Save as PDF
pdf('Top 10 CC Enrichment Dot Plot.pdf', width = 7, height = 5)

ggplot(top_significant_CC_major, aes(x = reorder(Description, -p.adjust), y = p.adjust, size = Count, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 Cellular Component (CC) Enrichment Dot Plot", x = "Cellular Component", y = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

dev.off()

# Save as PNG
png('Top 10 CC Enrichment Dot Plot.png', width = 7, height = 5, units = "in", res = 300)

ggplot(top_significant_CC_major, aes(x = reorder(Description, -p.adjust), y = p.adjust, size = Count, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 Cellular Component (CC) Enrichment Dot Plot", x = "Cellular Component", y = "Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())

dev.off()



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################




significant_GO_gwas <- data.frame(matrix(NA, nrow = 0, ncol = 14))
colnames(significant_GO_gwas) <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", 
                                    "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", 
                                    "geneID", "Count", "GO", 'Type') 

significant_KEGG_gwas <- data.frame(matrix(NA, nrow = 0, ncol = 15))
colnames(significant_KEGG_gwas) <- c("category", "subcategory", "ID", "Description", "GeneRatio", 
                                      "BgRatio", "RichFactor", "FoldEnrichment", "zScore", "pvalue", 
                                      "p.adjust", "qvalue", "geneID", "Count", "Type")





########### all protein coding loci in GWAS data

protein_coding <- data[data$Gene.Type %like% 'protein coding',] %>%
  separate_rows(Gene.Type, Gene.Name, Type, sep = ",") %>%
  .[.$Gene.Type %like% 'protein coding',]

protein_coding$Gene.Name <- gsub(" ", "", protein_coding$Gene.Name)

write.table(unique(protein_coding$Gene.Name),'all protein coding loci.txt',row.names = F,
            col.names = F, sep = '\t', quote = F)

gene_ids <- bitr(unique(protein_coding$Gene.Name), fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("protein_coding", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("protein_coding", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("protein_coding", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)

write.xlsx(GO,'all protein coding loci GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)


# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("protein_coding", nrow(kegg_dataframe))

write.xlsx(as.data.frame(kegg),'all protein coding loci KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

##############
##############
######### Score


score <- data.frame(mean_mutation = NA,
                    count_deg = NA,
                    group_size= NA,
                    group= NA)


score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% protein_coding$Gene.Name,]$total.adj),
                              length(intersect(protein_coding$Gene.Name,COAD_all$Gene)),
                              length(protein_coding$Gene.Name),
                              "protein_coding")


############### all protein coding loci up


all_up_protein <- data[data$Distance!=0 & data$Gene.Name==0,]

all_up_protein_gene_name <- gsub(' ','',all_up_protein$Maped.downe.gene) %>%
  unique()

write.table(all_up_protein_gene_name,'all up protein coding loci.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(all_up_protein_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("all_up_protein", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("all_up_protein", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("all_up_protein", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("all_up_protein", nrow(kegg_dataframe))


write.xlsx(as.data.frame(kegg),'all up protein coding loci KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)


######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% all_up_protein$Maped.downe.gene,]$total.adj),
                              length(intersect(all_up_protein$Maped.downe.gen,COAD_all$Gene)),
                              length(all_up_protein$Maped.downe.gen),
                              "all_up_protein")


############### all protein coding loci up to 2000

up_to_2000 <- all_up_protein[all_up_protein$Distance <= 2000,]
up_to_2000_gene_name <- gsub(' ','',up_to_2000$Maped.downe.gene) %>%
  unique()

write.table(up_to_2000_gene_name,'all up protein coding loci up to 2000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_2000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_2000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_2000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_2000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to 2000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_2000", nrow(kegg_dataframe))


write.xlsx(as.data.frame(kegg),'all up protein coding loci up to 2000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_2000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_2000$Maped.downe.gen,COAD_all$Gene)),
                              length(up_to_2000$Maped.downe.gen),
                              "up_to_2000")


############### all protein coding loci up to 10000

up_to_10000 <- all_up_protein[all_up_protein$Distance <= 10000,]
up_to_10000_gene_name <- gsub(' ','',up_to_10000$Maped.downe.gene) %>%
  unique()

write.table(up_to_10000_gene_name,'all up protein coding loci up to 10000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_10000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_10000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_10000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_10000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to 10000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_10000", nrow(kegg_dataframe))

write.xlsx(as.data.frame(kegg),'all up protein coding loci up to 10000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_10000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_10000$Maped.downe.gen,COAD_all$Gene)),
                              length(up_to_10000$Maped.downe.gen),
                              "up_to_10000")


############### all protein coding loci up to 20000

up_to_20000 <- all_up_protein[all_up_protein$Distance <= 20000,]
up_to_20000_gene_name <- gsub(' ','',up_to_20000$Maped.downe.gene) %>%
  unique()

write.table(up_to_20000_gene_name,'all up protein coding loci up to 20000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_20000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_20000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_20000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_20000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)

write.xlsx(GO,'all up protein coding loci up to 20000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_20000", nrow(kegg_dataframe))

write.xlsx(as.data.frame(kegg),'all up protein coding loci up to 20000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_20000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_20000$Maped.downe.gen,COAD_all$Gene)),
                              length(up_to_20000$Maped.downe.gen),
                              "up_to_20000")


############### all protein coding loci up to 50000

up_to_50000 <- all_up_protein[all_up_protein$Distance <= 50000,]
up_to_50000_gene_name <- gsub(' ','',up_to_50000$Maped.downe.gene) %>%
  unique()

write.table(up_to_50000_gene_name,'all up protein coding loci up to 50000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_50000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_50000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_50000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_50000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to 50000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_50000", nrow(kegg_dataframe))


write.xlsx(as.data.frame(kegg),'all up protein coding loci up to 50000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_50000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_50000$Maped.downe.gen, COAD_all$Gene)),
                              length(up_to_50000$Maped.downe.gen),
                              "up_to_50000")

############### all protein coding loci up to 100000

up_to_100000 <- all_up_protein[all_up_protein$Distance <= 100000,]
up_to_100000_gene_name <- gsub(' ','',up_to_100000$Maped.downe.gene) %>%
  unique()

write.table(up_to_100000_gene_name,'all up protein coding loci up to 100000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_100000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_100000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_100000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_100000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to 100000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_100000", nrow(kegg_dataframe))


write.xlsx(as.data.frame(kegg),'all up protein coding loci up to 100000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_100000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_100000$Maped.downe.gen, COAD_all$Gene)),
                              length(up_to_100000$Maped.downe.gen),
                              "up_to_100000")


######## all protein coding loci up to between 2000-10000

up_to_between_2000_10000 <- all_up_protein[all_up_protein$Distance > 2000 & all_up_protein$Distance <=10000,]
up_to_between_2000_10000_gene_name <- gsub(' ','',up_to_between_2000_10000$Maped.downe.gene) %>%
  unique()

write.table(up_to_between_2000_10000_gene_name,'all up protein coding loci up to between 2000-10000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_between_2000_10000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_between_2000_10000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_between_2000_10000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_between_2000_10000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to between 2000-10000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_between_2000_10000", nrow(kegg_dataframe))

write.xlsx(as.data.frame(kegg),'all up protein coding loci up to between 2000-10000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_between_2000_10000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_between_2000_10000$Maped.downe.gen, COAD_all$Gene)),
                              length(up_to_between_2000_10000$Maped.downe.gen),
                              "up_to_between_2000_10000")

######## all protein coding loci up to between 10000-20000

up_to_between_10000_20000 <- all_up_protein[all_up_protein$Distance > 10000 & all_up_protein$Distance <=20000,]
up_to_between_10000_20000_gene_name <- gsub(' ','',up_to_between_10000_20000$Maped.downe.gene) %>%
  unique()

write.table(up_to_between_10000_20000_gene_name,'all up protein coding loci up to between 10000-20000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_between_10000_20000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_between_10000_20000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_between_10000_20000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_between_10000_20000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to between 10000-20000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_between_10000_20000", nrow(kegg_dataframe))


write.xlsx(as.data.frame(kegg),'all up protein coding loci up to between 10000-20000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_between_10000_20000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_between_10000_20000$Maped.downe.gen,COAD_all$Gene)),
                              length(up_to_between_10000_20000$Maped.downe.gen),
                              "up_to_between_10000_20000")


######## all protein coding loci up to between 20000-50000

up_to_between_20000_50000 <- all_up_protein[all_up_protein$Distance > 20000 & all_up_protein$Distance <=50000,]
up_to_between_20000_50000_gene_name <- gsub(' ','',up_to_between_20000_50000$Maped.downe.gene) %>%
  unique()

write.table(up_to_between_20000_50000_gene_name,'all up protein coding loci up to between 20000-50000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_between_20000_50000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_between_20000_50000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_between_20000_50000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_between_20000_50000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to between 20000-50000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_between_20000_50000", nrow(kegg_dataframe))


write.xlsx(as.data.frame(kegg),'all up protein coding loci up to between 20000-50000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_between_20000_50000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_between_20000_50000$Maped.downe.gen, COAD_all$Gene)),
                              length(up_to_between_20000_50000$Maped.downe.gen),
                              "up_to_between_20000_50000")



######## all protein coding loci up to between 50000-100000

up_to_between_50000_100000 <- all_up_protein[all_up_protein$Distance > 50000 & all_up_protein$Distance <=100000,]
up_to_between_50000_100000_gene_name <- gsub(' ','',up_to_between_50000_100000$Maped.downe.gene) %>%
  unique()

write.table(up_to_between_50000_100000_gene_name,'all up protein coding loci up to between 50000-100000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_between_50000_100000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_between_50000_100000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_between_50000_100000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_between_50000_100000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to between 50000-100000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_between_50000_100000", nrow(kegg_dataframe))


write.xlsx(as.data.frame(kegg),'all up protein coding loci up to between 50000-100000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_between_50000_100000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_between_50000_100000$Maped.downe.gen,COAD_all$Gene)),
                              length(up_to_between_50000_100000$Maped.downe.gen),
                              "up_to_between_50000_100000")


######## all protein coding loci up to higher than 100000

up_to_higher_than_100000 <- all_up_protein[all_up_protein$Distance > 100000,]
up_to_higher_than_100000_gene_name <- gsub(' ','',up_to_higher_than_100000$Maped.downe.gene) %>%
  unique()

write.table(up_to_higher_than_100000_gene_name,'all up protein coding loci up to higher than 100000.txt',row.names = F,quote = F,
            col.names = F,sep = '\t')

gene_ids <- bitr(up_to_higher_than_100000_gene_name, fromType="SYMBOL", toType="ENTREZID", org.Hs.eg.db)

BP <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

MF <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

CC <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

BP_dataframe <- as.data.frame(BP)
BP_dataframe$GO <- rep("BP", nrow(BP_dataframe))
BP_dataframe$Type <- rep("up_to_higher_than_100000", nrow(BP_dataframe))


MF_dataframe <- as.data.frame(MF)
MF_dataframe$GO <- rep("MF", nrow(MF_dataframe))
MF_dataframe$Type <- rep("up_to_higher_than_100000", nrow(MF_dataframe))

CC_dataframe <- as.data.frame(CC)
CC_dataframe$GO <- rep("CC",nrow(CC_dataframe))
CC_dataframe$Type <- rep("up_to_higher_than_100000", nrow(CC_dataframe))

GO <- rbind(BP_dataframe,MF_dataframe)
GO <- rbind(GO,CC_dataframe)


write.xlsx(GO,'all up protein coding loci up to higher than 100000 GO.xlsx',row.names = F)
significant_GO_gwas <- rbind(significant_GO_gwas,GO)

# KEGG pathway enrichment analysis
kegg <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_dataframe <- as.data.frame(kegg)
kegg_dataframe$Type <- rep("up_to_higher_than_100000", nrow(kegg_dataframe))


write.xlsx(as.data.frame(kegg),'all up protein coding loci up to higher than 100000 KEGG.xlsx',row.names = F)
significant_KEGG_gwas <- rbind(significant_KEGG_gwas, kegg_dataframe)

######### Score

score[nrow(score) + 1, ] <- c(mean(crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% up_to_higher_than_100000$Maped.downe.gene,]$total.adj),
                              length(intersect(up_to_higher_than_100000$Maped.downe.gen,COAD_all$Gene)),
                              length(up_to_higher_than_100000$Maped.downe.gen),
                              "up_to_higher_than_100000")



final_score <- score
final_score  <- na.omit(final_score )
final_score$DEGs_percent <- round(as.numeric(final_score$count_deg)/as.numeric(final_score$group_size),4)*100
final_score$mean_mutation <- as.numeric(final_score$mean_mutation)%>%
  round(.,5)
final_score$count_deg <- as.numeric(final_score$count_deg)


final_score <- subset(final_score, select= c("group","group_size","count_deg",
                                             "DEGs_percent","mean_mutation"))


mutation_data <- list()
mutation_data[["protein_coding"]] <- crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% protein_coding$Gene.Name,]$total.adj
group_data <- final_score$group

for (group in group_data[-1]) {
  
  df <- get(group)

  mutation_data [[group]]<- crad_gene_mutation[crad_gene_mutation$Hugo_Symbol %in% df$Maped.downe.gene,]$total.adj
  }

# Create a dataframe from all groups
mutation_data_df <- data.frame(
  value = c(mutation_data$protein_coding, 
            mutation_data$all_up_protein,
            mutation_data$up_to_2000,
            mutation_data$up_to_10000,
            mutation_data$up_to_20000,
            mutation_data$up_to_50000,
            mutation_data$up_to_100000,
            mutation_data$up_to_between_2000_10000,
            mutation_data$up_to_between_10000_20000,
            mutation_data$up_to_between_20000_50000,
            mutation_data$up_to_between_50000_100000,
            mutation_data$up_to_higher_than_100000),
  
  group = rep(c("protein_coding", 
                "all_up_protein",
                "up_to_2000",
                "up_to_10000",
                "up_to_20000",
                "up_to_50000",
                "up_to_100000",
                "up_to_between_2000_10000",
                "up_to_between_10000_20000",
                "up_to_between_20000_50000",
                "up_to_between_50000_100000",
                "up_to_higher_than_100000"), 
              c(length(mutation_data$protein_coding), 
                length(mutation_data$all_up_protein),
                length(mutation_data$up_to_2000),
                length(mutation_data$up_to_10000),
                length(mutation_data$up_to_20000),
                length(mutation_data$up_to_50000),
                length(mutation_data$up_to_100000),
                length(mutation_data$up_to_between_2000_10000),
                length(mutation_data$up_to_between_10000_20000),
                length(mutation_data$up_to_between_20000_50000),
                length(mutation_data$up_to_between_50000_100000),
                length(mutation_data$up_to_higher_than_100000)))
)



# Perform ANOVA
anova_result <- aov(value ~ group, data = mutation_data_df)

# Print the ANOVA summary
summary(anova_result)

## PostHOC
tukey_result <- TukeyHSD(anova_result)
summary(tukey_result)
View(as.data.frame(tukey_result$group))

posthoc_bonferroni <- pairwise.t.test(mutation_data_df$value, mutation_data_df$group, 
                                      p.adjust.method = "bonferroni")

View(posthoc_bonferroni$p.value)


# Define a vector of 12 colors
colors <- c("skyblue", "lightgreen", "yellow", "pink", "lightblue", "orange", "purple", 
            "red", "blue", "green", "brown", "gray")

legend_names <- c("All Upstream", "Protein Coding", "0–10,000 bp", "0–100,000 bp", "0–2,000 bp", 
                  "0–20,000 bp", "0–50,000 bp", "10,000–20,000 bp", "2,000–10,000 bp", "20,000–50,000 bp", 
                  "50,000–100,000 bp", "More than 100,000 bp")

outlier_df <- data.frame(group = "protein_coding", value = 0.27, label = "*")

# Create the box plot
pdf('box plot.pdf')


ggplot(mutation_data_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = colors, labels = legend_names) + 
  theme_minimal() +
  labs(title = "Distribution of Adjusted Mutation Counts by Group", y = "Value", fill = "Legend") +
  ylim(0, 0.27) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) +
  geom_text(data = outlier_df, aes(x = group, y = value, label = label), 
            color = "red", size = 6, vjust = -0.5) +
  geom_text(data = outlier_df, aes(x = group, y = value, label = "1.01"), 
            color = "red", size = 4, vjust = 0.5) 

dev.off()



#### PNG
# Your plot code

plot <- ggplot(mutation_data_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = colors, labels = legend_names) + 
  theme_minimal() +
  labs(title = "Distribution of Adjusted Mutation Counts by Group", y = "Value", fill = "Legend") +
  ylim(0, 0.27) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) +
  geom_text(data = outlier_df, aes(x = group, y = value, label = label), 
            color = "red", size = 6, vjust = -0.5) +
  geom_text(data = outlier_df, aes(x = group, y = value, label = "1.01"), 
            color = "red", size = 4, vjust = 0.5)

# Save the plot as a PNG file
ggsave("box plot.png", plot = plot, width = 8, height = 6, dpi = 300)



############ Mutation Vs Distance

subset_data <- subset(data, select = c("SNPS", "Gene.Type", "Gene.Name", 
                                       "Maped.downe.gene", "Distance", "Type"))

subset_data <- subset_data %>%
  separate_rows(Gene.Type, Gene.Name, Type, sep = ",")

subset_data$Gene.Type <- gsub(' ','',subset_data$Gene.Type)
subset_data$Gene.Name <- gsub(' ','',subset_data$Gene.Name)
subset_data$Type <- gsub(' ','',subset_data$Type)
subset_data$Maped.downe.gene <- gsub(' ','',subset_data$Maped.downe.gene)
subset_data <- subset_data[subset_data$Gene.Type %in% c('proteincoding',0),]

### Plot

plot_mutation_distance <- data.frame(name = NA,
                                     distance = NA,
                                     mutation = NA)

# Loop through subset_data
for (i in 1:nrow(subset_data)) {
  
  # Check if Gene.Type is not 0
  if (subset_data$Gene.Type[i] != 0) {
    
    # Get mutation or NA if empty
    mutation_value <- crad_gene_mutation[crad_gene_mutation$Hugo_Symbol == subset_data$Gene.Name[i], ]$total.adj
    mutation_value <- ifelse(length(mutation_value) == 0, NA, mutation_value)
    
    # Add row
    plot_mutation_distance[nrow(plot_mutation_distance) + 1, ] <- c(
      subset_data$Gene.Name[i],
      0,
      mutation_value
    )
    
  } else {
    
    # Get mutation or NA if empty
    mutation_value <- crad_gene_mutation[crad_gene_mutation$Hugo_Symbol == subset_data$Maped.downe.gene[i], ]$total.adj
    mutation_value <- ifelse(length(mutation_value) == 0, NA, mutation_value)
    
    # Add row
    plot_mutation_distance[nrow(plot_mutation_distance) + 1, ] <- c(
      subset_data$Maped.downe.gene[i],
      subset_data$Distance[i],
      mutation_value
    )
  }
}

plot_mutation_distance <- plot_mutation_distance[-1,]
plot_mutation_distance <- plot_mutation_distance[!duplicated(plot_mutation_distance),]
plot_mutation_distance <- na.omit(plot_mutation_distance)
plot_mutation_distance$distance <- as.numeric(plot_mutation_distance$distance)
plot_mutation_distance$mutation <- as.numeric(plot_mutation_distance$mutation)

plot_mutation_distance <- plot_mutation_distance %>%
  group_by(name) %>%
  slice_min(distance) %>%
  ungroup()

# Linear Model
lm_model <- lm(mutation ~ distance, data = plot_mutation_distance)
summary(lm_model)  # Check p-value, R², etc.

### Draw Plot

pdf('mutation Vs distance.pdf')

ggplot(data = plot_mutation_distance, aes(x = distance, y = mutation)) +
  geom_point(color = "blue", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linetype = "dashed") +  # Regression Line
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x = 0.5, label.y = 0.27, size = 3) +  # Equation, R², p-value
  labs(title = "Distance vs Mutation with Regression", 
       x = "Distance", 
       y = "Mutation Count") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(plot_mutation_distance$distance, na.rm = TRUE),
                                  max(plot_mutation_distance$distance, na.rm = TRUE), 
                                  length.out = 7),
                     labels = number_format(accuracy = 1)) +
  scale_y_continuous(breaks = seq(0, 0.284, length.out = 7),   
                     labels = number_format(accuracy = 0.001),
                     limits = c(0, 0.284)) +
  annotate("text", x = 0, y = 0.284, label = "*", size = 6, color = "red") +  # Star for outlier
  annotate("text", x = 0, y = 0.280, label = "1.014", size = 3, color = "red") + # Outlier label
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16))

dev.off()

######### Chi-Square

chi_square_data <- subset(final_score, select = c("group", "group_size", "count_deg"))
chi_square_data$group_size <- as.numeric(chi_square_data$group_size)
chi_square_data$count_deg <- as.numeric(chi_square_data$count_deg)
chi_square_data$group_size <- chi_square_data$group_size - chi_square_data$count_deg
colnames(chi_square_data) <- c('Groups','non DEGs','DEGs')
rownames(chi_square_data) <-chi_square_data$Groups
chi_square_data <- chi_square_data[,-1]


# Perform the Chi-Square Test
chi_result <- chisq.test(chi_square_data)

# Print the results
print(chi_result)





#### Hypergeometric test

degs_hyper_results_polymorphism_groups <- data.frame(group= NA,
                                                     group_gene= NA,
                                                     DEGs= NA,
                                                     hyper_p_value= NA)


for (i in c(1:nrow(final_score))){

  hyper_p_value <- phyper(q = final_score$count_deg[i] - 1,      
                          m = as.numeric(final_score$group_size[i]),           
                          n = 25000 - as.numeric(final_score$group_size[i]), 
                          k = length(COAD_all$Gene),          
                          lower.tail = FALSE) 
  
  degs_hyper_results_polymorphism_groups[nrow(degs_hyper_results_polymorphism_groups)+1,] <- c(final_score$group[i],
                                                                                               final_score$group_size[i],
                                                                                               final_score$count_deg[i],
                                                                                               hyper_p_value)
}












######################
######################
######################


#### Create Manhattan plot for all loci

data_for_manhattan <- read.xlsx('colorectal cancer final GWAS edited.xlsx',sheetIndex = 1)

data_for_manhattan <- subset(data_for_manhattan, select= c("SNPS","Gene.Type","Gene.Name","Maped.downe.gene",
                                                           "Distance","Type","CHR_ID","CHR_POS","REGION","P.VALUE"))

colnames(data_for_manhattan) <- c("SNP","Gene.Type","Gene.Name","Maped.downe.gene",
                                  "Distance","Type","CHR","BP","REGION","P")


# Calculate -log10(P-value)
data_for_manhattan$P <- as.numeric(data_for_manhattan$P)
data_for_manhattan$logP <- -log10(data_for_manhattan$P)

### Numeric BP
data_for_manhattan[data_for_manhattan$SNP=='rs71167281',]$BP <- '222034660'
data_for_manhattan[data_for_manhattan$SNP=='rs67052019',]$BP <- '109822839'
data_for_manhattan$BP <- as.numeric(data_for_manhattan$BP)


# Order the data by CHR and BP
data_for_manhattan <- data_for_manhattan[order(data_for_manhattan$CHR, data_for_manhattan$BP), ]


### Get chromosome sizes for GRCh38 (hg38)
chr_lengths <- SeqinfoForUCSCGenome("hg38") %>%
  as.data.frame()

chr_lengths <- chr_lengths[c(1:23),]

chr_lengths <- data.frame(
  CHR = rownames(chr_lengths),
  chr_length = chr_lengths$seqlengths)

chr_lengths$CHR <- gsub('chr','',chr_lengths$CHR)

# Calculate cumulative start positions based on chromosome lengths
chr_lengths$cumulative_start <- c(0, cumsum(as.numeric(chr_lengths$chr_length[-length(chr_lengths$chr_length)])))


# Merge the chromosome lengths with the GWAS data
data_for_manhattan <- data_for_manhattan %>%
  left_join(chr_lengths, by = "CHR") %>%
  mutate(CHR_POS = BP + cumulative_start)



# Convert CHR to factor for proper ordering
data_for_manhattan$CHR <- factor(data_for_manhattan$CHR, levels = c(1:22, "X"))


# Calculate the midpoint for each chromosome for labeling the x-axis
chr_ticks <- chr_lengths %>%
  mutate(chr_center = cumulative_start + chr_length / 2)



######## Plotting the Manhattan Plot
pdf('manhattan plot.pdf', width = 15, height = 10)

ggplot(data_for_manhattan, aes(x = CHR_POS, y = logP, color = as.factor(CHR))) +
  geom_point(alpha = 0.6, size = 1) +  # Points for each SNP
  scale_color_manual(values = rep(c("blue", "red"), length.out = length(unique(data_for_manhattan$CHR)))) +  # Alternating blue and red colors
  scale_x_continuous(
    breaks = chr_ticks$chr_center,  # Position for the x-axis tick labels
    labels = chr_ticks$CHR  # Chromosome labels
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Angling the x-axis labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16)  # Center the title and increase its size
  ) +
  labs(
    x = "Chromosome", 
    y = "-log10(P-value)",
    title = "Manhattan Plot with Proportional Chromosome Scaling",
    color = "Chromosome"
  )



dev.off()





######## Plotting the Manhattan Plot with outliers
pdf('manhattan_plot_with_outliers.pdf', width = 15, height = 10)

# Identify outliers (y > 65)
outliers <- data_for_manhattan %>% filter(logP > 65)

# Adjust the main data to limit y-axis
data_for_manhattan_limited <- data_for_manhattan %>%
  mutate(logP = ifelse(logP > 65, 65, logP))  # Cap y-axis values at 60

# Plot the Manhattan plot
ggplot(data_for_manhattan_limited, aes(x = CHR_POS, y = logP, color = as.factor(CHR))) +
  geom_point(alpha = 0.6, size = 1) +  # Points for each SNP
  geom_point(data = outliers, aes(x = CHR_POS, y = 65), shape = 8, size = 3, color = "black") +  # Add stars for outliers
  scale_color_manual(values = rep(c("blue", "red"), length.out = length(unique(data_for_manhattan$CHR)))) +  # Alternating blue and red colors
  scale_x_continuous(
    breaks = chr_ticks$chr_center,  # Position for the x-axis tick labels
    labels = chr_ticks$CHR  # Chromosome labels
  ) +
  scale_y_continuous(
    limits = c(0, 65),  # Limit y-axis to 0-60
    expand = expansion(mult = c(0, 0.05))  # Add space above for stars
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Angling the x-axis labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16)  # Center the title and increase its size
  ) +
  labs(
    x = "Chromosome", 
    y = "-log10(P-value)",
    title = "Manhattan Plot with Outliers Highlighted",
    color = "Chromosome"
  )

dev.off()





#### Create Manhattan plot for protein coding loci

data_for_manhattan <- read.xlsx('colorectal cancer final GWAS edited.xlsx',sheetIndex = 1)

data_for_manhattan <- data_for_manhattan[data_for_manhattan$Gene.Type %like% 'protein coding',]

data_for_manhattan <- subset(data_for_manhattan, select= c("SNPS","Gene.Type","Gene.Name","Maped.downe.gene",
                                                           "Distance","Type","CHR_ID","CHR_POS","REGION","P.VALUE"))

colnames(data_for_manhattan) <- c("SNP","Gene.Type","Gene.Name","Maped.downe.gene",
                                  "Distance","Type","CHR","BP","REGION","P")


# Calculate -log10(P-value)
data_for_manhattan$P <- as.numeric(data_for_manhattan$P)
data_for_manhattan$logP <- -log10(data_for_manhattan$P)

### Numeric BP
#data_for_manhattan[data_for_manhattan$SNP=='rs71167281',]$BP <- '222034660'
#data_for_manhattan[data_for_manhattan$SNP=='rs67052019',]$BP <- '109822839'
data_for_manhattan$BP <- as.numeric(data_for_manhattan$BP)


# Order the data by CHR and BP
data_for_manhattan <- data_for_manhattan[order(data_for_manhattan$CHR, data_for_manhattan$BP), ]


### Get chromosome sizes for GRCh38 (hg38)
chr_lengths <- SeqinfoForUCSCGenome("hg38") %>%
  as.data.frame()

chr_lengths <- chr_lengths[c(1:23),]

chr_lengths <- data.frame(
  CHR = rownames(chr_lengths),
  chr_length = chr_lengths$seqlengths)

chr_lengths$CHR <- gsub('chr','',chr_lengths$CHR)

# Calculate cumulative start positions based on chromosome lengths
chr_lengths$cumulative_start <- c(0, cumsum(as.numeric(chr_lengths$chr_length[-length(chr_lengths$chr_length)])))


# Merge the chromosome lengths with the GWAS data
data_for_manhattan <- data_for_manhattan %>%
  left_join(chr_lengths, by = "CHR") %>%
  mutate(CHR_POS = BP + cumulative_start)



# Convert CHR to factor for proper ordering
data_for_manhattan$CHR <- factor(data_for_manhattan$CHR, levels = c(1:22, "X"))


# Calculate the midpoint for each chromosome for labeling the x-axis
chr_ticks <- chr_lengths %>%
  mutate(chr_center = cumulative_start + chr_length / 2)



######## Plotting the Manhattan Plot
pdf('manhattan plot protein coding.pdf', width = 15, height = 10)

ggplot(data_for_manhattan, aes(x = CHR_POS, y = logP, color = as.factor(CHR))) +
  geom_point(alpha = 0.6, size = 1) +  # Points for each SNP
  scale_color_manual(values = rep(c("blue", "red"), length.out = length(unique(data_for_manhattan$CHR)))) +  # Alternating blue and red colors
  scale_x_continuous(
    breaks = chr_ticks$chr_center,  # Position for the x-axis tick labels
    labels = chr_ticks$CHR  # Chromosome labels
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Angling the x-axis labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16)  # Center the title and increase its size
  ) +
  labs(
    x = "Chromosome", 
    y = "-log10(P-value)",
    title = "Manhattan Plot with Proportional Chromosome Scaling",
    color = "Chromosome"
  )



dev.off()




#### Create Manhattan plot for non coding loci

data_for_manhattan <- read.xlsx('colorectal cancer final GWAS edited.xlsx',sheetIndex = 1)

data_for_manhattan <- data_for_manhattan[data_for_manhattan$Maped.downe.gene!=0,]

data_for_manhattan <- subset(data_for_manhattan, select= c("SNPS","Gene.Type","Gene.Name","Maped.downe.gene",
                                                           "Distance","Type","CHR_ID","CHR_POS","REGION","P.VALUE"))

colnames(data_for_manhattan) <- c("SNP","Gene.Type","Gene.Name","Maped.downe.gene",
                                  "Distance","Type","CHR","BP","REGION","P")


# Calculate -log10(P-value)
data_for_manhattan$P <- as.numeric(data_for_manhattan$P)
data_for_manhattan$logP <- -log10(data_for_manhattan$P)

### Numeric BP
data_for_manhattan[data_for_manhattan$SNP=='rs71167281',]$BP <- '222034660'
data_for_manhattan[data_for_manhattan$SNP=='rs67052019',]$BP <- '109822839'
data_for_manhattan$BP <- as.numeric(data_for_manhattan$BP)


# Order the data by CHR and BP
data_for_manhattan <- data_for_manhattan[order(data_for_manhattan$CHR, data_for_manhattan$BP), ]


### Get chromosome sizes for GRCh38 (hg38)
chr_lengths <- SeqinfoForUCSCGenome("hg38") %>%
  as.data.frame()

chr_lengths <- chr_lengths[c(1:23),]

chr_lengths <- data.frame(
  CHR = rownames(chr_lengths),
  chr_length = chr_lengths$seqlengths)

chr_lengths$CHR <- gsub('chr','',chr_lengths$CHR)

# Calculate cumulative start positions based on chromosome lengths
chr_lengths$cumulative_start <- c(0, cumsum(as.numeric(chr_lengths$chr_length[-length(chr_lengths$chr_length)])))


# Merge the chromosome lengths with the GWAS data
data_for_manhattan <- data_for_manhattan %>%
  left_join(chr_lengths, by = "CHR") %>%
  mutate(CHR_POS = BP + cumulative_start)



# Convert CHR to factor for proper ordering
data_for_manhattan$CHR <- factor(data_for_manhattan$CHR, levels = c(1:22, "X"))


# Calculate the midpoint for each chromosome for labeling the x-axis
chr_ticks <- chr_lengths %>%
  mutate(chr_center = cumulative_start + chr_length / 2)



######## Plotting the Manhattan Plot
pdf('manhattan plot non coding.pdf', width = 15, height = 10)

ggplot(data_for_manhattan, aes(x = CHR_POS, y = logP, color = as.factor(CHR))) +
  geom_point(alpha = 0.6, size = 1) +  # Points for each SNP
  scale_color_manual(values = rep(c("blue", "red"), length.out = length(unique(data_for_manhattan$CHR)))) +  # Alternating blue and red colors
  scale_x_continuous(
    breaks = chr_ticks$chr_center,  # Position for the x-axis tick labels
    labels = chr_ticks$CHR  # Chromosome labels
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Angling the x-axis labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16)  # Center the title and increase its size
  ) +
  labs(
    x = "Chromosome", 
    y = "-log10(P-value)",
    title = "Manhattan Plot with Proportional Chromosome Scaling",
    color = "Chromosome"
  )



dev.off()




write.xlsx(significant_KEGG_major,'significant KEGG in major band.xlsx',row.names = F)
write.xlsx(significant_GO_major,'significant GO in major band.xlsx',row.names = F)
write.xlsx(mutual_significant_major_band_KEGG,'significant KEGG mutual.xlsx',row.names = F)
write.xlsx(mutual_significant_major_band_GO,'significant GO mutual.xlsx',row.names = F)
write.xlsx(significant_KEGG_major[significant_KEGG_major$ID %in% c('hsa05033', 'hsa04742', 'hsa04740'),],
           'significant KEGG survival.xlsx',row.names = F)
write.xlsx(significant_GO_major[significant_GO_major$ID %in% c("GO:0001580", "GO:0050909", "GO:0001906", "GO:0033038", 
                                                               "GO:0008527", "GO:0140375", "GO:0070821", "GO:0070820", 
                                                               "GO:0043086", "GO:0051346", "GO:0019730","GO:0004857",
                                                               "GO:1902894", "GO:0070085", "GO:0097014", "GO:0036126"),],
           'significant GO survival.xlsx',row.names = F)








rs_data_20q13 <- data.frame("Gencode.Id" = NA,
                            "Gene.Symbol" = NA,
                            "Variant.Id" = NA,
                            "SNP.Id" = NA,
                            "P.Value" = NA,
                            "NES" = NA,
                            "Tissue" = NA)
for (i in c(1:22)){
  
  temporary <- read.csv(paste0('New folder/GTEx Portal (',i ,')','.csv'))
  rs_data_20q13 <- rbind(rs_data_20q13,temporary)
}


rs_data_20q13 <- na.omit(rs_data_20q13)


ensembl_ids <- c("ENSG00000232043", "ENSG00000270951", "ENSG00000225806", 
                 "ENSG00000273619", "ENSG00000275437", "ENSG00000226332", 
                 "ENSG00000233017")


gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = ensembl)


unique(rs_data_20q13$Gene.Symbol)


gene_20q_13 <- read.table('20q13 genes.txt')


intersect(gene_20q_13$V1,COAD_all$Gene)





# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic",
               "Non-Coding", "Non-coding (annotated)", "Non-coding (uncharacterized)"),
  Count = c(10, 1, 1, 9, 1,8),
  Group = c("All", "Coding", "Coding-sub", 
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[4,'Percentage'] <- round(plot_data[4,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "Non-coding (annotated)", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_10p14.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_10p14.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()


###################################
###################################
###################################
###################################




# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic",
               "Non-Coding", "Non-coding (annotated)", "Non-coding (uncharacterized)"),
  Count = c(18, 16, 16, 2, 1,1),
  Group = c("All", "Coding", "Coding-sub", 
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[4,'Percentage'] <- round(plot_data[4,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "Non-coding (annotated)", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_10q25.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_10q25.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()


###################################
###################################
###################################
###################################



# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", '3 prime UTR',
               "Non-Coding", "Non-coding (uncharacterized)"),
  Count = c(13, 6, 5, 1, 7, 7),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[5,'Percentage'] <- round(plot_data[5,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "3 prime UTR", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "3 prime UTR" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_11q12.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_11q12.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()



###################################
###################################
###################################
###################################


# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", "Splice Acceptor",
               "Non-Coding",'Non-coding (annotated)', "Non-coding (uncharacterized)"),
  Count = c(21, 11, 10, 1, 10, 6,4),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[5,'Percentage'] <- round(plot_data[5,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "Splice Acceptor", 'Non-coding (annotated)',"Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Splice Acceptor" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_12p13.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_12p13.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()



###################################
###################################
###################################
###################################



# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic",
               "Non-Coding","Non-coding (uncharacterized)"),
  Count = c(18, 10, 10, 8, 8),
  Group = c("All", "Coding", "Coding-sub", 
            "Non-Coding", "Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[4,'Percentage'] <- round(plot_data[4,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding",
                                        "Intronic", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_15q13.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_15q13.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()


###################################
###################################
###################################
###################################



# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", "3 prime UTR",
               "Non-Coding",'Non-coding (annotated)', "Non-coding (uncharacterized)"),
  Count = c(20, 15, 14, 1, 5, 2,3),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[5,'Percentage'] <- round(plot_data[5,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "3 prime UTR", 'Non-coding (annotated)',"Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Splice Acceptor" = "#D633FF",  # Electric Purple for Others
            "3 prime UTR" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_18q21.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_18q21.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()




###################################
###################################
###################################
###################################



# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", "Missense","Synonymous",
               "Non-Coding",'Non-coding (annotated)', "Non-coding (uncharacterized)"),
  Count = c(26, 18, 16, 1,1, 8, 2,6),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[6,'Percentage'] <- round(plot_data[6,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "Missense","Synonymous", 'Non-coding (annotated)',
                                        "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Synonymous" = "#D633FF",  # Electric Purple for Others
            "3 prime UTR" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_19q13.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_19q13.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()




###################################
###################################
###################################
###################################




# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic",
               "Non-Coding", "Non-coding (annotated)", "Non-coding (uncharacterized)"),
  Count = c(18, 1, 1, 17, 10,7),
  Group = c("All", "Coding", "Coding-sub", 
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[4,'Percentage'] <- round(plot_data[4,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "Non-coding (annotated)", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_1q41.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_1q41.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()


###################################
###################################
###################################
###################################



# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic",
               "Non-Coding", "Non-coding (annotated)", "Non-coding (uncharacterized)"),
  Count = c(27, 3, 3, 24, 2,22),
  Group = c("All", "Coding", "Coding-sub", 
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[4,'Percentage'] <- round(plot_data[4,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "Non-coding (annotated)", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_20p12.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_20p12.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()


###################################
###################################
###################################
###################################




# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", "3 prime UTR",
               "Non-Coding",'Non-coding (annotated)', "Non-coding (uncharacterized)"),
  Count = c(29, 16, 14, 2, 13, 3,10),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[5,'Percentage'] <- round(plot_data[5,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "3 prime UTR", 'Non-coding (annotated)',"Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Splice Acceptor" = "#D633FF",  # Electric Purple for Others
            "3 prime UTR" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_20q13.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_20q13.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()




###################################
###################################
###################################
###################################



# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic",
               "Non-Coding", "Non-coding (annotated)", "Non-coding (uncharacterized)"),
  Count = c(27, 9, 9, 18, 3,15),
  Group = c("All", "Coding", "Coding-sub", 
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[4,'Percentage'] <- round(plot_data[4,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "Non-coding (annotated)", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_6p21.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_6p21.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()


###################################
###################################
###################################
###################################





# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", "5 prime UTR",
               "Non-Coding",'Non-coding (annotated)', "Non-coding (uncharacterized)"),
  Count = c(32, 6, 5, 1, 26, 17,9),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[5,'Percentage'] <- round(plot_data[5,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "5 prime UTR", 'Non-coding (annotated)',"Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Splice Acceptor" = "#D633FF",  # Electric Purple for Others
            "5 prime UTR" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_8q24.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_8q24.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()




###################################
###################################
###################################
###################################




# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", "3 prime UTR",
               "Non-Coding",'Non-coding (annotated)', "Non-coding (uncharacterized)"),
  Count = c(18, 8, 5, 3, 10, 2,8),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub","Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[5,'Percentage'] <- round(plot_data[5,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic",
                                        "3 prime UTR", 'Non-coding (annotated)',"Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Splice Acceptor" = "#D633FF",  # Electric Purple for Others
            "3 prime UTR" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_9q34.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_9q34.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()




###################################
###################################
###################################
###################################





c('10p14, 10q25, 11q12, 12p13, 15q13, 18q21, 19q13, 1q41, 20p12, 20q13, 6p21, 8q24, 9q34') 




table(data[data$major.band=='9q34',]$Gene.Type)
table(data[data$major.band=='9q34',]$Type)






# Create a data frame with hierarchical structure
plot_data <- data.frame(
  Category = c("Polymorphisms", "Coding", "Intronic", "Missense", 'UTRs (5 prime UTR and 3 prime UTR)', "Others (Synonymous, Splice Acceptor, and Stop Gained)",
               "Non-Coding", "Non-coding (annotated)", "Non-coding (uncharacterized)"),
  Count = c(1346, 632, 577, 20, 24, 11, 714, 211, 503),
  Group = c("All", "Coding", "Coding-sub", "Coding-sub", "Coding-sub", "Coding-sub",
            "Non-Coding", "Non-Coding-sub", "Non-Coding-sub")  # Defining the parent group
)

# Compute percentage
plot_data <- plot_data %>% 
  group_by(Group) %>% 
  mutate(Percentage = round((Count / sum(Count[Group != "All"])) * 100, 1))

# Assign ring levels (Polymorphisms = Inner Circle, Coding and Non-Coding = Middle Circle, others = Outer Circle)
plot_data <- plot_data %>% 
  mutate(Ring = case_when(
    Category == "Polymorphisms" ~ "Inner",
    Group %in% c("Coding", "Non-Coding") ~ "Middle",
    TRUE ~ "Outer"
  ))

plot_data[1,'Percentage'] <- 100.0
plot_data[2,'Percentage'] <- round(plot_data[2,'Count'] / plot_data[1,'Count']*100, 1)
plot_data[7,'Percentage'] <- round(plot_data[7,'Count'] / plot_data[1,'Count']*100, 1)

# Reorder the categories manually for correct plotting order
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Polymorphisms", "Coding", "Non-Coding", "Intronic", "Missense", 
                                        "Others (Synonymous, Splice Acceptor, and Stop Gained)", 
                                        "UTRs (5 prime UTR and 3 prime UTR)", "Non-coding (annotated)", "Non-coding (uncharacterized)"))

# Define colors for each category, with red shades for Coding and blue shades for Non-Coding
colors <- c("Polymorphisms" = "#FFFFFF",           # White for Polymorphisms (Center)
            "Coding" = "#5E0075",                 # Deep Purple for Coding
            "Intronic" = "#8B00A3",               # Vivid Violet for Intronic
            "Missense" = "#B300E0",               # Strong Magenta for Missense
            "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",  # Electric Purple for Others
            "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",                   # Soft Pinkish Purple for UTRs
            "Non-Coding" = "#005F73",             # Deep Teal for Non-Coding
            "Non-coding (annotated)" = "#008C99", # Ocean Blue-Teal for Non-coding (annotated)
            "Non-coding (uncharacterized)" = "#00C4CC"  # Bright Cyan for Non-coding (uncharacterized)
)

# Save as PDF
pdf("nested_pie_chart_three_levels_updated_color_v5.pdf")

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))  # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PDF device
dev.off()

# Save as PNG with high resolution
png("nested_pie_chart_three_levels_updated_color_v5_high_quality.png", 
    width = 3000, height = 3000, res = 300)  # Higher resolution (300 DPI) for better quality

# Create Nested Pie Chart (Donut Chart) WITHOUT percentages
ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = NA) +  # Remove gaps by setting color = NA
  coord_polar(theta = "y", start = 0) +  # Convert to circular form
  scale_fill_manual(values = colors) +  # Use the updated color palette
  theme_void() +  # Remove background and axis
  theme(legend.title = element_blank(),
        legend.key = element_blank(),  # Remove legend keys
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title
        plot.margin = margin(0, 0, 0, 0))   # Remove outer margin
#+ labs(title = "Genomic Distribution of GWAS CRC-Associated Polymorphisms")  # New title

# Close the PNG device
dev.off()


