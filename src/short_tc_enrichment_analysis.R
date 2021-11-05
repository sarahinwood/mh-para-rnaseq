library('data.table')
library('fgsea')
library('ggplot2')

trinotate_report <- fread("data/mh_edited_transcript_ids/trinotate_annotation_report.txt", na.strings = ".")
gene_ids <- trinotate_report[!is.na(gene_ontology_pfam), unique(`#gene_id`)]
res_group <- fread("output/mh_timecourse/short_tc_deseq2/m120_res_group.csv")

go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_pfam, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))

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
fwrite(sorted_fgsea_res, "output/mh_timecourse/short_tc_fgsea/fgsea_developing_wasp_120vsC.csv")

##subset into only sig terms and merge w/annotations
sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.05)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/mh_timecourse/short_tc_fgsea/sig_annot_fgsea_pfam.csv")
##split into 3 tables --> biological process, cellular component and molecular function
bp_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="biological_process"]
cc_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="cellular_component"]
mf_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="molecular_function"]

##plot normalised enrichment for GO terms where padj<0.1 (but indicate if padj<0.05) - can change to only bp, cc or mf
ggplot(bp_res, aes(reorder(pathway_name, NES), NES)) +
  geom_text(aes(label=round(padj, digits=3)), vjust=0, hjust=0) +
  geom_col() +
  coord_flip() +
  labs(x="Biological Process GO Pathway", y="FGSEA Normalized Enrichment Score") + 
  theme_minimal()+
  theme(axis.text.y = element_text(size=20), axis.title = element_text(size=15))

####Core members that contribute to ES score SIGNAL TRANSDUCTION (present in list before running sum reaches max.dev. from 0)
sig_trans_res <- fgsea_res[fgsea_res$pathway == "GO:0007165",]
sig_trans_leading_edge <- data.frame(sig_trans_res$leadingEdge)
setnames(sig_trans_leading_edge, old=c("c..TRINITY_DN3907_c0_g1....TRINITY_DN453_c1_g4....TRINITY_DN1729_c0_g1..."), new=c("gene_id"))
sig_trans_leading_annots <- merge(sig_trans_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(sig_trans_leading_annots, "output/short_tc_fgsea/sig_trans/sig_trans_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0007165"]], ranks) + labs(title="signal transduction")

##CORE MEMBERS OF TRANSCRIPTION REGULATION
trans_reg_res <- fgsea_res[fgsea_res$pathway == "GO:0006355",]
trans_reg_leading_edge <- data.frame(trans_reg_res$leadingEdge)
setnames(trans_reg_leading_edge, old=c("c..TRINITY_DN2808_c0_g1....TRINITY_DN466_c4_g1....TRINITY_DN358_c9_g1..."), new=c("gene_id"))
trans_reg_leading_annots <- merge(trans_reg_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(trans_reg_leading_annots, "output/short_tc_fgsea/trans_reg/trans_reg_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0006355"]], ranks) + labs(title="regulation of transcription, DNA templated")

##CORE MEMBERS OF BIOSYNTHETIC
biosynth_res <- fgsea_res[fgsea_res$pathway == "GO:0009058",]
biosynth_leading_edge <- data.frame(biosynth_res$leadingEdge)
setnames(biosynth_leading_edge, old=c("c..TRINITY_DN14084_c0_g1....TRINITY_DN2014_c0_g1....TRINITY_DN292_c0_g1..."), new=c("gene_id"))
biosynth_leading_annots <- merge(biosynth_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(biosynth_leading_annots, "output/short_tc_fgsea/biosynth/biosynth_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0009058"]], ranks) + labs(title="biosynthetic process")
