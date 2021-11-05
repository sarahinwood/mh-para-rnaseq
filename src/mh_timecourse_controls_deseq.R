library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("EnhancedVolcano")
library("VennDiagram")

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

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)

######remove samples to see how many genes are still DE######
##remove 8h NC2&3 - find sample number in table
colData(dds)
##remove 8hNC2 - 4th in table
dds <- dds[,-4]
##remove 8hNC3 - sample 5 once 8hNC2 removed
dds <- dds[,-5]

#####OR####remove 8hNC1
dds <- dds[,-2]

##check correct samples were removed
colData(dds)
#############################################################

##create dds object for group analysis
dds_group <- copy(dds)
##create groupings of tissue+treatment
dds_group$group <- factor(paste(dds$Tissue,dds$Treatment,dds$Flow_cell,sep="_"))
##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)

resultsNames(dds_group)

##Make table of results for mh_timecourse_controls vs control heads
res_group <- results(dds_group, contrast = c("group", "Abdomen_Control_CDT11ANXX", "Abdomen_Control_CCU1EANXX"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/mh_timecourse_controls/deseq2/res_group.csv")
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/mh_timecourse_controls/deseq2/mh_timecourse_controls_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)
##Sub in any gene of interest to plot counts  
plotCounts(dds_group, "MH_TRINITY_DN13_c1_g1", intgroup = c("group"), main="")
##volcano plot
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)

##read in annotated transcriptome
trinotate_report <- fread("data/mh_edited_transcript_ids/trinotate_annotation_report.txt")
setnames(ordered_sig_res_group_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
sig_w_annots <- merge (ordered_sig_res_group_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/mh_timecourse_controls/deseq2/sig_genes_with_annots.csv")


only_8h1 <- fread("output/mh_timecourse_controls/deseq2/8hNC1_only_mh_timecourse_controls_analysis_sig_degs.csv")
only_8h1_degs <- only_8h1$rn

wout_8h1 <- fread("output/mh_timecourse_controls/deseq2/8hNC2&3_only_mh_timecourse_controls_analysis_sig_degs.csv")
wout_8h1_degs <- wout_8h1$rn

vd <- venn.diagram(x = list("8hNC1 vs old NCs"=only_8h1_degs, "8hNC2&3 vs old NCs"=wout_8h1_degs), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

##compare old vs new NC DEGs to para DEGs
old_vs_new_NCs <- fread("output/mh_timecourse_controls/deseq2/mh_timecourse_controls_analysis_sig_degs.csv")
old_vs_new_NCs$`#gene_id` <- tstrsplit(old_vs_new_NCs$rn, "MH_", keep=c(2))
old_vs_new_NC_ids <- old_vs_new_NCs$`#gene_id`
mh_para_m240 <- fread("short_tc_output/mh_timecourse/deseq2/control_vs_m240_sig_degs.csv")
mh_para_m240_ids <- mh_para_m240$rn

vd <- venn.diagram(x = list("Old vs New NC DEGs"=old_vs_new_NC_ids, "m240 DEGs"=mh_para_m240_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)
