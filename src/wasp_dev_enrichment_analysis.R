library('data.table')
library('fgsea')
library('ggplot2')

trinotate_report <- fread("data/trinotate_annotation_report.txt", na.strings = ".")
gene_ids <- trinotate_report[!is.na(gene_ontology_pfam), unique(`#gene_id`)]
res_group <- fread("output/mh_timecourse/deseq2/control_vs_m120_all_degs.csv")

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

go_term_list <- lapply(gene_ids, EXTRACT_GO_TERMS, trinotate_report=trinotate_report)
go_term_table <- rbindlist(go_term_list)
term_to_gene <- go_term_table[,list(list(gene_id)), by=accessions]
pathways <- term_to_gene[,V1]
names(pathways) <- term_to_gene[,accessions]

##use stat column from deseq results to rank genes (can change if wanted)
setorder(res_group, stat)
ranks <- res_group[!is.na(stat), stat]
names(ranks) <- res_group[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
fwrite(sorted_fgsea_res, "output/fgsea/fgsea_developing_wasp_120vsC.csv")

##read in file with functions added to GO terms when padj<0.1
annot_fgsea_res <- fread("output/fgsea/annot_fgsea_developing_wasp_120vsC.csv")
##split into 3 tables --> biological process, cellular component and molecular function
bp_res <- annot_fgsea_res[annot_fgsea_res$`pathway_kind`=="biological process"]
cc_res <- annot_fgsea_res[annot_fgsea_res$`pathway_kind`=="cellular component"]
mf_res <- annot_fgsea_res[annot_fgsea_res$`pathway_kind`=="molecular function"]

##plot normalised enrichment for GO terms where padj<0.1 (but indicate if padj<0.05) - can change to only bp, cc or mf
ggplot(cc_res, aes(reorder(pathway_name, NES), NES)) +
  geom_text(aes(label=round(padj, digits=3)), vjust=0, hjust=0) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Cellular Component GO Pathway", y="FGSEA Normalized Enrichment Score") + 
  theme_minimal()

