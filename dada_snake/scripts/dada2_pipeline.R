library(dada2)
library(ggplot2)
library(dplyr)
library(tidyr)
# Snakemake passes params, input, output objects automatically
input_dir <- snakemake@input[[1]]
truncLen  <- as.integer(snakemake@params[["truncLen"]])
trimLeft  <- as.integer(snakemake@params[["trimLeft"]])
min_len   <- as.integer(snakemake@params[["min_len"]])
max_len   <- as.integer(snakemake@params[["max_len"]])

#rds objects and quality control stats
out_seqtab  <- snakemake@output[["seqtab"]]
out_qcplot  <- snakemake@output[["qc_plot"]]
out_qcplot_filtered <- snakemake@output[["qc_filter_plot"]]
out_qcstats <- snakemake@output[["qc_stats"]]
# plots 
out_reads_plot <- snakemake@output[["reads_plot"]]
out_errF_plot  <- snakemake@output[["errF_plot"]]
out_errR_plot  <- snakemake@output[["errR_plot"]]

log_file <- snakemake@log[[1]]
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

progress <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
  cat(msg, file = stderr())   # → visible in Snakemake terminal
  message(msg)                # → also captured in log file
}

# --- File discovery ---
fnFs <- sort(list.files(input_dir, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(input_dir, pattern = "_2.fastq", full.names = TRUE))
sample.names <- sub("_1.fastq", "", basename(fnFs))

# Guard: stop early with a clear message if no files found
if (length(fnFs) == 0) {
  stop("No forward FASTQ files found in: ", input_dir,
       "\nCheck input_dir path and file naming pattern (_1.fastq.gz)")
}

message("Found ", length(fnFs), " samples: ", paste(sample.names, collapse = ", "))

# --- Quality plots ---
pdf(out_qcplot)
plotQualityProfile(fnFs[1:min(4, length(fnFs))])
plotQualityProfile(fnRs[1:min(4, length(fnRs))])
dev.off()

# --- Filter & trim ---
filtFs <- file.path(input_dir, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(input_dir, "filtered", paste0(sample.names, "_R_filt.fastq"))

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen  = truncLen,
  trimLeft  = trimLeft,
  maxN = 0, maxEE = c(2, 2),
  truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = TRUE
)

pdf(out_qcplot_filtered)
plotQualityProfile(filtFs[file.exists(filtFs)][1:min(4, length(filtFs[file.exists(filtFs)]))])
plotQualityProfile(filtRs[file.exists(filtRs)][1:min(4, length(filtRs[file.exists(filtRs)]))])
dev.off()

# --- QC stats ---
qc_stats <- data.frame(
  sample   = sample.names,
  reads_in = out[, 1],
  reads_out = out[, 2],
  fraction_passing = out[, 2] / out[, 1]
)
write.table(qc_stats, out_qcstats, sep = "\t", row.names = FALSE, quote = FALSE)

pdf(out_reads_plot, width = 10, height = 5)
barplot_data <- t(as.matrix(out[, c("reads.in", "reads.out")]))
colnames(barplot_data) <- sample.names
barplot(
  barplot_data,
  beside    = TRUE,
  col       = c("#4E9A8A", "#C0392B"),
  las       = 2,
  cex.names = 0.7,
  main      = "Read counts before and after filtering",
  ylab      = "Number of reads",
  legend    = c("Before filtering", "After filtering"),
  args.legend = list(x = "topright", bty = "n")
)
dev.off()
progress("Filtering complete — writing QC stats")
# --- Learn errors ---
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
progress("Error models learned (errF + errR)")

p_errF <- plotErrors(errF, nominalQ = TRUE)
ggsave(out_errF_plot, plot = p_errF, width = 10, height = 8)

p_errR <- plotErrors(errR, nominalQ = TRUE)
ggsave(out_errR_plot, plot = p_errR, width = 10, height = 8)
progress("Plots saved")
# --- Denoise ---
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
progress("Denoising complete")
# --- Merge ---
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
progress("Mergers complete")

# --- Seqtab ---
seqtab <- makeSequenceTable(mergers)
progress("Seqtab ready")

# --- Remove chimeras ---
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
progress("Removed chimeras")

# --- Subset by amplicon length ---
seq_lengths <- nchar(colnames(seqtab.nochim))
seqtab.nochim <- seqtab.nochim[, seq_lengths >= min_len & seq_lengths <= max_len]
progress("Subsetting complete")

saveRDS(seqtab.nochim, out_seqtab)
progress("Seqtab RDS saved")

if (ncol(seqtab.nochim) == 0) {
  stop(
    "seqtab.nochim has 0 ASVs after length subsetting (", min_len, "-", max_len, " bp).\n",
    "ASV lengths present before subsetting: ",
    paste(sort(unique(nchar(colnames(seqtab)))), collapse = ", "), "\n",
    "Widen min_length/max_length in your YAML config."
  )
}

cat("[INFO] ", ncol(seqtab.nochim), " ASVs passed length filter\n", file = stderr())
    
taxa <- assignTaxonomy(seqtab.nochim, snakemake@params[["silva_db"]], multithread = TRUE)
saveRDS(taxa, snakemake@output[["taxa"]])
progress("Taxonomy assigned and RDS saved")

sink(type = "message")
sink(type = "output")
close(log_con)