library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)
library(VennDiagram)
library(svglite)

gene2tx <- fread("data/mh_edited_transcript_ids/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files from salmon filtering because res. from filtering with STAR and lost bro that way
quant_files <- list.files(path="output/asw_mh_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_key.csv")
setkey(sample_data, Sample_name)

trinotate_report <- fread("data/mh_edited_transcript_ids/trinotate_annotation_report.txt", na.strings=".")

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
counts_matrix <- counts(dds)

##look at counts for pcr gene
dds_parasitism_PCR <- copy(dds)
##create groupings of weevil location and parasitism response to parasitoid
dds_parasitism_PCR$group <- factor(paste(dds$Parasitism_PCR))
##plot counts for PCR target gene
plotCounts(dds_parasitism_PCR, "MH_TRINITY_DN5286_c0_g1", intgroup = c("group"))

######remove samples to see how many genes are still DE######
##remove 8h samples
colData(dds)
##repeat 6x to remove all 8h sample at top of table
dds <- dds[,-1]
dds <- dds[,-1]
dds <- dds[,-1]
dds <- dds[,-1]
dds <- dds[,-1]
dds <- dds[,-1]
##check correct samples were removed
colData(dds)
#############################################################

##create dds object for group analysis
dds_group <- copy(dds)
##create groupings of tissue+treatment
dds_group$group <- factor(paste(dds$Tissue,dds$Treatment,sep="_"))
##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
saveRDS(dds_group, file = "output/mh_timecourse/short_tc_deseq2/dds_group.rds")

dds_group <- readRDS("output/mh_timecourse/short_tc_deseq2/dds_group.rds")
resultsNames(dds_group)

##Make table of results - m30
m30_res_group <- results(dds_group, contrast = c("group", "Abdomen_m30", "Abdomen_Control"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
m30_ordered_res_group <- m30_res_group[order(m30_res_group$padj),]
##Make data table and write to output
m30_ordered_res_group_table <- data.table(data.frame(m30_ordered_res_group), keep.rownames = TRUE)
m30_ordered_sig_res_group_table <- subset(m30_ordered_res_group_table, padj < 0.05)
m30_sig_annots <- merge(m30_ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id")
fwrite(m30_ordered_res_group_table, "output/mh_timecourse/short_tc_deseq2/m30_res_group.csv")
fwrite(m30_ordered_sig_res_group_table, "output/mh_timecourse/short_tc_deseq2/m30_mh_timecourse_controls_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)
fwrite(m30_sig_annots, "output/mh_timecourse/short_tc_deseq2/m30_sig_annots.csv")

##Make table of results - m120
m120_res_group <- results(dds_group, contrast = c("group", "Abdomen_m120", "Abdomen_Control"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
m120_ordered_res_group <- m120_res_group[order(m120_res_group$padj),]
##Make data table and write to output
m120_ordered_res_group_table <- data.table(data.frame(m120_ordered_res_group), keep.rownames = TRUE)
m120_ordered_sig_res_group_table <- subset(m120_ordered_res_group_table, padj < 0.05)
m120_sig_annots <- merge(m120_ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id")
fwrite(m120_ordered_res_group_table, "output/mh_timecourse/short_tc_deseq2/m120_res_group.csv")
fwrite(m120_ordered_sig_res_group_table, "output/mh_timecourse/short_tc_deseq2/m120_mh_timecourse_controls_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)
fwrite(m120_sig_annots, "output/mh_timecourse/short_tc_deseq2/m120_sig_annots.csv")

##Make table of results - m240
m240_res_group <- results(dds_group, contrast = c("group", "Abdomen_m240", "Abdomen_Control"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
m240_ordered_res_group <- m240_res_group[order(m240_res_group$padj),]
##Make data table and write to output
m240_ordered_res_group_table <- data.table(data.frame(m240_ordered_res_group), keep.rownames = TRUE)
m240_ordered_sig_res_group_table <- subset(m240_ordered_res_group_table, padj < 0.05)
m240_sig_annots <- merge(m240_ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id")
fwrite(m240_ordered_res_group_table, "output/mh_timecourse/short_tc_deseq2/m240_res_group.csv")
fwrite(m240_ordered_sig_res_group_table, "output/mh_timecourse/short_tc_deseq2/m240_mh_timecourse_controls_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)
fwrite(m240_sig_annots, "output/mh_timecourse/short_tc_deseq2/m240_sig_annots.csv")

##Sub in any gene of interest to plot counts  
plot_gene <- plotCounts(dds_group, "MH_TRINITY_DN7765_c1_g1", 
                        intgroup = c("Sample_name"), returnData = TRUE)
ggplot(plot_gene,
       aes(x = Sample_name, y = count)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() + xlab("Parasitism Timepoint") + ylab("Normalized Count")

##volcano plot
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)

m30_names <- m30_ordered_sig_res_group_table$rn
fwrite(list(m30_names), "m30_names.txt")
m120_names <- m120_ordered_sig_res_group_table$rn
fwrite(list(m120_names), "m120_names.txt")
m240_names <- m240_ordered_sig_res_group_table$rn
fwrite(list(m240_names), "m240_names.txt")

Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("30 minute DEGs"=m30_names, "120 minute DEGs"=m120_names, "240 minute DEGs"=m240_names),
                   filename=NULL, fill=Set1, alpha=0.4, cex = 1, cat.cex=1, cat.dist = 0.05,
                   cat.pos = 180, lwd=1)
grid.newpage()
grid.draw(vd)
##save venn diagram svg so can edit in inkscape
ggsave(file="output/mh_timecourse/short_tc_deseq2/short_tc_venn.svg", plot=vd, width=10, height=4)

##what genes aren't sig at 240m
m120_not_m240 <- setdiff(m120_names, m240_names)
m30_not_m240 <- setdiff(m30_names, m240_names)
not_sig_m240 <- data.table(unique(c(m120_not_m240, m30_not_m240)))
not_sig_m240_deseq_res <- merge(not_sig_m240, m240_ordered_res_group_table, by.x="V1", by.y="rn")
not_sig_m240_annots <- merge(not_sig_m240_deseq_res, trinotate_report, by.x="V1", by.y="#gene_id")
fwrite(not_sig_m240_annots, "output/mh_timecourse/short_tc_deseq2/m240_not_sig.csv")




