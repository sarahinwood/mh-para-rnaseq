library("tximport")
library("data.table")
library("DESeq2")

##ASW gene to trans map (edited to have ID's match concat transcriptome - ASW_TRINITY)
gene2tx <- fread("data/mh_edited_transcript_ids/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files - quantified against concat ASW-MH transcriptome
quant_files <- list.files(path="output/asw_mh_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_key.csv", header=TRUE)
setkey(sample_data, Sample_name)


##Create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
saveRDS(dds, file = "output/mh_timecourse/deseq2/dds.rds")
dds <- readRDS("output/mh_timecourse/deseq2/dds.rds")

dds$parasitism_pcr <- factor(paste(dds$Parasitism_PCR))
##plot expression pattern for gene
counts <- plotCounts(dds, "MH_TRINITY_DN1053_c1_g1", intgroup = c("parasitism_pcr"), main = "M. hyperodae Multiplex PCR Target", returnData=TRUE)

PCR_gene_counts <- plotCounts(dds, "MH_TRINITY_DN1053_c1_g1", intgroup = c("parasitism_pcr"), main = "M. hyperodae Multiplex PCR Target", returnData = TRUE)


##plot expression pattern for viral genes
plotCounts(dds, "MH_TRINITY_DN17931_c0_g1", intgroup = c("parasitism_pcr"), main = "Bro-N Domain containing gene #2")


##8h18ha has 15 count for PCR target

mh_longest_trinotate <- fread("data/mh_edited_transcript_ids/trinotate_longest_isoform.csv")
