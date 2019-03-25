library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("EnhancedVolcano")

gene2tx <- fread("data/mh_transcriptome/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

  ##Find all salmon quant files
quant_files <- list.files(path="output/mh_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
  ##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
  ##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
  ##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
  ##Import table describing samples
sample_data <- fread("data/sample_key.csv")
setkey(sample_data, Sample_name)

  ##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
  ##save dds object
saveRDS(dds, file = "output/exposed/deseq2/dds.rds")

 ##create dds object for group analysis
dds_group <- copy(dds)
  ##create groupings of tissue+treatment
dds_group$group <- factor(paste(dds$Tissue,dds$Treatment,sep="_"))
  ##add group to design
design(dds_group) <- ~group
  ##run deseq2 and generate results
dds_group <- DESeq(dds_group)
  ##save dds_group
saveRDS(dds_group, file = "output/exposed/deseq2/dds_group.rds")

resultsNames(dds_group)

  ##Make table of results for exposed vs control heads
res_group <- results(dds_group, contrast = c("group", "Head_Exposed", "Head_Control"), lfcThreshold = 1, alpha = 0.1)
  ##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
  ##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/exposed/deseq2/res_group.csv")
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/exposed/deseq2/exposed_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)
##Sub in any gene of interest to plot counts  
plotCounts(dds_group, "TRINITY_DN38122_c0_g1", intgroup = c("group"), main="bro [Spodoptera frugiperda granulovirus] (padj 1, L2FC 5.3)")
##volcano plot
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)

  ##read in annotated transcriptome
trinotate_report <- fread("data/mh_transcriptome/trinotate_annotation_report.txt")
setnames(ordered_sig_res_group_table, old=c("rn"), new=c("#gene_id"))
  ##merge list of sig genes with annotations
sig_w_annots <- merge (ordered_sig_res_group_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
  ##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/exposed/deseq2/sig_genes_with_annots.csv")
