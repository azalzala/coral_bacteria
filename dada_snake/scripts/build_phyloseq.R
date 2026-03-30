library(phyloseq)

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ── Load inputs ────────────────────────────────────────────
asv_count <- readRDS(snakemake@input[["seqtab"]])
taxa      <- readRDS(snakemake@input[["taxa"]])
metadata  <- read.csv(snakemake@params[["metadata"]], check.names = FALSE)

message("Dimensions on load:")
message("  asv_count: ", paste(dim(asv_count), collapse = " x "))
message("  taxa:      ", paste(dim(taxa),      collapse = " x "))
message("  metadata:  ", paste(dim(metadata),  collapse = " x "))

# ── Fix rownames ───────────────────────────────────────────
# Taxa rownames must match asv_count column names (ASV sequences)
rownames(taxa) <- colnames(asv_count)

# Metadata rownames must be sample IDs matching asv_count rownames
rownames(metadata) <- rownames(asv_count)

# ── Dimension sanity checks ────────────────────────────────
stopifnot(
  "Taxa rows must equal ASV count columns" =
    ncol(asv_count) == nrow(taxa),
  "Metadata rows must equal ASV count rows" =
    nrow(asv_count) == nrow(metadata)
)

# ── Reorder metadata to match asv_count sample order ──────
metadata <- metadata[rownames(asv_count), , drop = FALSE]

# ── Build phyloseq components ──────────────────────────────
otu_ps <- otu_table(asv_count, taxa_are_rows = FALSE)
tax_ps <- tax_table(as.matrix(taxa))
smd_ps <- sample_data(metadata)

# ── Name alignment checks ──────────────────────────────────
message("taxa_names match (otu vs tax): ",
        identical(taxa_names(otu_ps), taxa_names(tax_ps)))

message("sample_names match (otu vs smd): ",
        identical(sample_names(otu_ps), sample_names(smd_ps)))

if (!identical(taxa_names(otu_ps), taxa_names(tax_ps))) {
  stop("Taxa names do not match between OTU table and tax table. ",
       "Check that rownames(taxa) <- colnames(asv_count) was applied correctly.")
}

if (!identical(sample_names(otu_ps), sample_names(smd_ps))) {
  stop("Sample names do not match between OTU table and sample data. ",
       "Check that rownames(metadata) <- metadata$run_SRR matches rownames(asv_count).")
}

# ── Build phyloseq object ──────────────────────────────────
ps <- phyloseq(otu_ps, tax_ps, smd_ps)

message("Phyloseq object built successfully:")
print(ps)

saveRDS(ps, snakemake@output[["ps"]])

sink(type = "message")
sink(type = "output")
close(log_con)