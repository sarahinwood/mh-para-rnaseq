library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
#library("Biostrings")
library("dplyr")
library("VennDiagram")
library("EnhancedVolcano")

####SET-UP
gene2tx <- fread("data/mh_edited_transcript_ids/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

  ##Find all salmon quant files
quant_files <- list.files(path="output/asw_mh_concat_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
  ##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
  ##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
  ##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
  ##Import table describing samples
sample_data <- fread("data/sample_key.csv")
setkey(sample_data, Sample_name)

  ##Create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)

dds_concat_group_mh <- copy(dds)
##plot counts of COI gene
dds_concat_group_mh$group <- factor(paste(dds$Parasitism_PCR))
plotCounts(dds_concat_group_mh, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"))

##make counts matrix
counts_matrix <- counts(dds)
counts_colsum <- data.table(data.frame(colSums(counts_matrix)), keep.rownames = TRUE)
fwrite(counts_colsum, "output/counts_colsum.csv")

####FILTER for genes with some counts
  ##Keep only genes with reasonable no. of reads (allows for more statistical power in later analysis)
kept_genes<-rownames(counts_matrix[rowMeans(counts_matrix)>5 | rowMax(counts_matrix)>10,])
  ##use to subset into KEPT GENES for later on in pairwise analysis
dds_group<-dds_group[kept_genes,]

####LRT analysis
##Select only abdomen samples
dds_abdo <- dds[kept_genes,dds$Tissue == "Abdomen"&dds$Treatment != "Control"]
  ##convert to factors
dds_abdo$Treatment <- factor(dds_abdo$Treatment)
dds_abdo$Wasp_Location <- factor(dds_abdo$Wasp_Location)
  ##add factors of ineterst to design
design(dds_abdo) <- ~Wasp_Location+Treatment
  ##run deseq, must specify reduced model for LRT test
dds_abdo <- DESeq(dds_abdo, test = "LRT", reduced = ~Wasp_Location)
  ##filter results on alpha
dds_abdo_res <- results(dds_abdo, alpha = 0.1)
 ##order based off padj
ordered_dds_abdo_res <- dds_abdo_res[order(dds_abdo_res$padj),]
  ##make list of sig genes
timecourse_sig_genes <- subset(dds_abdo_res, pvalue < 0.1)
  ##make list of sig gene names
sig_gene_names <- row.names(timecourse_sig_genes)
  ##save list of sig gene names
fwrite(data.table(sig_gene_names), "output/mh_timecourse/deseq2/timecourse_sig_gene_names.csv")

  ##Order results based of padj
ordered_sig_degs <- timecourse_sig_genes[order(timecourse_sig_genes$padj),]
  ##make datatable and write to output
timecourse_ord_degs_table <- data.table(data.frame(ordered_sig_degs), keep.rownames = TRUE)
fwrite(timecourse_ord_degs_table, "output/mh_timecourse/deseq2/timecourse_analysis_sig_degs.csv")
setnames(ordered_degs_table, old=c("rn"), new=c("#gene_id"))

##merge list of sig genes with annotations
sig_w_annots <- merge(ordered_degs_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/mh_timecourse/deseq2/sig_genes_with_annots.csv")


####PAIRWISE analysis
##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
saveRDS(dds, file = "output/mh_timecourse/deseq2/dds.rds")

##create dds object for group analysis
dds_group <- copy(dds)
##create groupings of tissue+treatment
dds_group$group <- factor(paste(dds$Tissue,dds$Treatment,sep="_"))
##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
##save dds_group
saveRDS(dds_group, file = "output/mh_timecourse/deseq2/dds_group.rds")

resultsNames(dds_group)
##Make table of results for exposed vs control heads
res_group <- results(dds_group, contrast = c("group", "Abdomen_m120", "Abdomen_Control"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
m120_ordered_res <- res_group[order(res_group$padj),]
##Make data table and write to output
m120_ordered_res_table <- data.table(data.frame(m120_ordered_res), keep.rownames = TRUE)
##filter for sig DEGs
m120_ordered_sig_res_table <- subset(m120_ordered_res_table, padj < 0.05)
fwrite(m120_ordered_sig_res_table, "output/deseq2/control_vs_m120_sig_degs.csv", col.names = TRUE, row.names = FALSE)

##volcano plot
EnhancedVolcano(m120_ordered_res_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)

##m120 vs control ALL - use for GO term enrichment analysis
fwrite(m120_ordered_res_table, "output/deseq2/control_vs_m120_all.csv", col.names = TRUE, row.names = FALSE)

m30_names <- m30_ordered_sig_res_table$rn
m120_names <- m120_ordered_sig_res_table$rn
m240_names <- m240_ordered_sig_res_table$rn

Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("30 minutes"=m30_names, "120 minutes"=m120_names, "240 minutes"=m240_names), filename=NULL, fill=Set1, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

unique(c(m30_names, m120_names, m240_names))

trinotate_report <- fread("data/mh_transcriptome/trinotate_annotation_report.txt", na.strings=".")
##merge list of sig genes with annotations
sig_w_annots <- merge(m120_ordered_sig_res_table, trinotate_report, by.x="rn", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/deseq2/control_v_120m_sig_genes_with_annots.csv")

##Compare DEG Annots from filtered and non-filtered
F_annots <- fread("output/deseq2/dedup_control_v_120m_sig_genes_with_annots.csv", na.strings = ".")
NF_annots <- fread("non_filtered_output/mh_timecourse/deseq2/dedup_control_v_120m_sig_genes_with_annots.csv", na.strings = ".")

F_annot_names <- F_annots[,tstrsplit(sprot_Top_BLASTX_hit, "^", fixed=TRUE, keep = 1)]
F_annot_noNA <- data.table(na.omit(F_annot_names$V1))
NF_annot_names <- NF_annots[,tstrsplit(sprot_Top_BLASTX_hit, "^", fixed=TRUE, keep = 1)]
NF_annot_noNA <- data.table(na.omit(NF_annot_names$V1))

vd <- venn.diagram(x = list("DEG Names (F)"=F_annot_noNA$V1, "DEG Names (NF)"=NF_annot_noNA$V1), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

unique_to_NF <- setdiff(NF_annot_noNA, F_annot_noNA)
fwrite(unique_to_NF, "output/NF_v_F/annots_unique_to_NF.csv")
unique_to_F <- setdiff(F_annot_noNA, NF_annot_noNA)
fwrite(unique_to_F, "output/NF_v_F/annots_unique_to_F.csv")
