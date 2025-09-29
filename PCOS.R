#This project involves differential gene analysis of Granulosa cells (PCOS vs Control), and annotation using GSEA

gc_counts <-read.csv("GSE155489_gc_pcos_counts.csv",row.names=1,check.names=FALSE)
sample <- c("GC_B7","GC_B8","GC_B15","GC_B16","GC_B13","GC_B14","GC_B2","GC_B30")
condition <- c("PCOS","PCOS","PCOS","PCOS","control","control","control","control")
gc_meta <- data.frame(sample,condition,stringsAsFactors = FALSE)
library(DESeq2)
ds_gc <-DESeqDataSetFromMatrix(
  countData= gc_counts,
  colData= gc_meta,
  design= ~ condition
)
deseq <- DESeq(ds_gc)
res_gc <- results(deseq,contrast = c("condition","PCOS","control"))
res_df <- as.data.frame(res_gc)
res_df <- na.omit(res_df)
gene_list <- res_df$log2FoldChange
names(gene_list) <- rownames(res_df)
gene_list <- sort(gene_list,decreasing = TRUE)
library(clusterProfiler) 
library(org.Hs.eg.db)
gsea_entrez <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list_entrez <- gene_list[names(gene_list) %in% gsea_entrez$SYMBOL]
names(gene_list_entrez) <- gsea_entrez$ENTREZID[match(names(gene_list_entrez),gsea_entrez$SYMBOL)]
gsea_go <- gseGO(geneList = gene_list_entrez,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",          # can be "BP", "MF", or "CC"
                 minGSSize = 10,
                 pvalueCutoff = 0.05,
                 verbose = TRUE)
gsea_go_df <- as.data.frame(gsea_go)

gsea_kegg <- gseKEGG(geneList = gene_list_entrez,
                     organism = "hsa",   # hsa = human
                     minGSSize = 10,
                     pvalueCutoff = 0.05,
                     verbose = TRUE)
gsea_kegg_df <- as.data.frame(gsea_kegg)