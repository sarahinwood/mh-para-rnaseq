#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)

###########
# GLOBALS #
###########

unmapped_names <- snakemake@input[["unmapped_names"]]

########
# MAIN #
########

unmapped_names <- read.delim(unmapped_names, header = FALSE)
fixed_names <- tstrsplit(unmapped_names$V1, " ", keep = 1)
fwrite(fixed_names, quote=FALSE, snakemake@output[["fixed_names"]])

#write log
sessionInfo()