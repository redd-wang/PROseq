#read in feature count output by Subread, followed by grep to consolidate read count including zero
CR0 <- read.table("35CR_counts_0min.txt", header = TRUE, row.names = 1)
CR60 <- read.table("35CR_counts_60min.txt", header = TRUE, row.names = 1)
CR90 <- read.table("35CR_counts_90min.txt", header = TRUE, row.names = 1)
CR120 <- read.table("35CR_counts_120min.txt", header = TRUE, row.names = 1)

#extract column 6, the gene count, to one count table
Lu35CR <- data.frame(CR0[ ,6], CR60[ ,6], CR90[ ,6], CR120[ ,6], row.names = rownames(CR0))
colnames(Lu35CR) <- c("X0min","X60min","X90min", "X120min")

#library load
library(org.Hs.eg.db)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(pheatmap)
library(AnnotationDbi)
library(biomaRt)

#construct colData and condition for DESeq2
col_df <- data.frame(colnames(Lu35CR))
condition <- factor(c("baseline", "60minposttx", "90minposttx", "120minposttx"))
col_df$condition <- condition

##the following is for no replicate data

#normalize count number
dds_normal_Lu35CR <- DESeqDataSetFromMatrix(countData = Lu35CR, colData = col_df, design = ~1)
dds_normal_Lu35CR <- estimateSizeFactors(dds_normal_Lu35CR)
Lu35CR_normalized <- counts(dds_normal_Lu35CR, normalize=TRUE)
Lu35CR_normalized <- as.data.frame(Lu35CR_normalized)

#cleanup Ensembl ID and align with gene symbol
Ensemblidwithversion <- row.names(Lu35CR)
Ensemblidwithoutversion <- sub("\\..*","",Ensemblidwithversion)
rownames(Lu35CR)<- Ensemblidwithoutversion
rownames(Lu35CR_normalized)<- Ensemblidwithoutversion
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),filters = "ensembl_gene_id",
  values = Ensemblidwithoutversion,
  mart = ensembl)
gene_map_unique <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
gene_map_ordered <- gene_map_unique[match(Ensemblidwithoutversion, gene_map_unique$ensembl_gene_id), ]
Lu35CR_normalized$genemap <- gene_map_ordered$hgnc_symbol

#calculate log2FC between 60min and 0min
pseudocount=1
Lu35CR_normalized$log2FC_60min <- log2((Lu35CR_normalized$X60min+pseudocount)/(Lu35CR_normalized$X0min +pseudocount))
rownames(Lu35CR_normalized) <- Ensemblidwithoutversion
Lu35CR_60minsorted <- Lu35CR_normalized[order(Lu35CR_normalized$log2FC_60min, decreasing = TRUE), ]

