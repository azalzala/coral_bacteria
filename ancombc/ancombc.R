library(ANCOMBC)
library(dplyr)
library(ggplot2)
library(purrr)
library(phyloseq)
library(tibble)
library(decontam)
library(metafor)
#renv::status()
#renv::restore()
#renv::snapshot()

ps1 <- readRDS('ps1.rds')
ps2 <- readRDS('ps2.rds')
ps3 <- readRDS('ps3.rds')
ps4 <- readRDS('ps4.rds')
ps5 <- readRDS('ps5.rds')
sample_data(ps1)$run_SRR #PRJNA341929
sample_data(ps2)$run_SRR # PRJNA967137
sample_data(ps3)$run_SRR #PRJNA510614
sample_data(ps4)$run_SRR # PRJNA1211024
sample_data(ps5)$run_SRR   # PRJNA1050834
sample_data(ps2)

# pre-processing phyloseq for merging
any(taxa_sums(ps2_order) ==0)
sum(taxa_sums(ps2_order) == 0)


# expand to genus in abundance measure
ps1_order <- tax_glom(ps1, taxrank = "Genus")
ps2_order <- tax_glom(ps2, taxrank = "Genus")
ps3_order <- tax_glom(ps3, taxrank = "Genus")
ps4_order <- tax_glom(ps4, taxrank = "Genus")
ps5_order <- tax_glom(ps5, taxrank = "Genus")


# remove taxons with less than one appearance across samples - unsuccessful for some reason 
ps1_clean<-prune_taxa(taxa_sums(ps1_order) > 0, ps1_order)
ps2_clean <-prune_taxa(taxa_sums(ps2_order) > 0, ps2_order)
ps3_clean<-prune_taxa(taxa_sums(ps3_order) > 0, ps3_order)
ps4_clean<-prune_taxa(taxa_sums(ps4_order) > 0, ps4_order)
ps5_clean<-prune_taxa(taxa_sums(ps5_order) > 0, ps5_order)

# Select only bacteria
ps1_order_ <- subset_taxa(ps1_clean, Kingdom == "Bacteria")
ps2_order_ <- subset_taxa(ps2_clean, Kingdom == "Bacteria")
ps3_order_ <- subset_taxa(ps3_clean, Kingdom == "Bacteria")
ps4_order_ <- subset_taxa(ps4_clean, Kingdom == "Bacteria")
ps5_order_ <- subset_taxa(ps5_clean, Kingdom == "Bacteria")

# Remove samples with less than 1000 reads

ps1_order_bacteria <- prune_samples(names(which(sample_sums(ps1_order_) >= 1000)), 
                    ps1_order_)
ps2_order_bacteria <- prune_samples(names(which(sample_sums(ps2_order_) >= 1000)), 
                                  ps2_order_)
ps3_order_bacteria <- prune_samples(names(which(sample_sums(ps3_order_) >= 1000)), 
                                  ps3_order_)
ps4_order_bacteria <- prune_samples(names(which(sample_sums(ps4_order_) >= 1000)), 
                                  ps4_order_)
ps5_order_bacteria <- prune_samples(names(which(sample_sums(ps5_order_) >= 1000)), 
                                  ps5_order_)


# Define your filtering steps as named phyloseq objects in order
ps_list <- list(
  "on import" = ps1,
  "Post: genus"    = ps1_order,
  "Post: taxon filter"          = ps1_clean,
  "Post: bacteria filter " = ps1_order_, 
  "Post: sample filter"      = ps1_order_bacteria
)


# summary table 
tracking <- tibble(
  Step        = names(ps_list),
  Taxa        = sapply(ps_list, ntaxa),
  Samples     = sapply(ps_list, nsamples)
) %>%
  mutate(
    Taxa_lost    = lag(Taxa) - Taxa,
    Pct_lost     = round(100 * Taxa_lost / lag(Taxa), 1)
  )

print(tracking)

### adding site columns for ps3
sample_data(ps3_order_bacteria)$site <- NA
grep("-[R]", sample_data(ps3_order_bacteria)$sample_id)
sample_data(ps3_order_bacteria)$site[c( 1,  3,  4, 14, 15, 16, 17, 18, 19, 20, 21, 22)] <- "Mangrove Lagoon"
sample_data(ps3_order_bacteria)$site[c(2, 5,  6,  7, 8,  9, 10, 11, 12, 13)] <- "Reef"


sample_data(ps5_order_bacteria)
# Mark seawater samples as negative controls in your metadata
sample_data(ps1_order_bacteria)$is_control <- sample_data(ps1_order_bacteria)$Species == "Seawater" 
sample_data(ps2_order_bacteria)$is_control <- sample_data(ps2_order_bacteria)$Species == "SW-21"  
sample_data(ps4_order_bacteria)$is_control <- sample_data(ps4_order_bacteria)$Species == "seawater"  
sample_data(ps5_order_bacteria)$is_control <- sample_data(ps5_order_bacteria)$Species == "seawater"

sample_data(ps1_order_bacteria)
sample_data(ps2_order_bacteria)
sample_data(ps3_order_bacteria)
sample_data(ps4_order_bacteria)
sample_data(ps5_order_bacteria)


# Prevalence method (recommended for low-biomass/eDNA)
contam_prev <- isContaminant(ps1_order_bacteria, method = "prevalence", neg = "is_control", threshold = 0.5)
contam_prev_2 <- isContaminant(ps2_order_bacteria, method = "prevalence", neg = "is_control", threshold = 0.5)
contam_prev_4 <- isContaminant(ps4_order_bacteria, method="prevalence", neg="is_control", threshold =0.5)
contam_prev_5 <- isContaminant(ps5_order_bacteria, method="prevalence", neg="is_control", threshold =0.5)

# See how many contaminants flagged
table(contam_prev$contaminant) # TRUE (116 / 116+373)
table(contam_prev_2$contaminant) # 0 
table(contam_prev_4$contaminant) # true 113/113+192 (leaves 192)
table(contam_prev_5$contaminant) # 172 true, 249 left (false)

ps1_order_bacteria <- prune_taxa(!contam_prev$contaminant, ps1_order_bacteria)
ntaxa(ps1_order_bacteria) # 373 taxa after filtering using seawater 

ps2_order_bacteria <- prune_taxa(!contam_prev_2$contaminant, ps2_order_bacteria)
ntaxa(ps2_order_bacteria) # 281

ps4_order_bacteria <- prune_taxa(!contam_prev_4$contaminant, ps4_order_bacteria)
ps5_order_bacteria <- prune_taxa(!contam_prev_5$contaminant, ps5_order_bacteria)


ps1_order_bacteria <- subset_samples(ps1_order_bacteria, Species != "Seawater")  
ps2_order_bacteria <- subset_samples(ps2_order_bacteria, Species != c("seawater")) 
ps2_order_bacteria <- subset_samples(ps2_order_bacteria, Species != c("negative_control")) 
ps2_order_bacteria <- subset_samples(ps2_order_bacteria, Species != c("artemia_salina")) 

ps4_order_bacteria <- subset_samples(ps4_order_bacteria, Species != "seawater")
ps4_order_bacteria <- subset_samples(ps4_order_bacteria, Species != "seawater")
ps4_order_bacteria <- subset_samples(ps4_order_bacteria, Species != "not_applicable")
ps5_order_bacteria <- subset_samples(ps5_order_bacteria, Species != c("seawater"))
sample_data(ps2_order_bacteria)$Species

# re-factor levels of the condition column so that LFC estimates are conditionHeat 
sample_data(ps1_order_bacteria)$condition <- factor(sample_data(ps1_order_bacteria)$condition, 
                                                    levels = c("normal", "heat"))
sample_data(ps2_order_bacteria)$condition <- factor(sample_data(ps2_order_bacteria)$condition, 
                                                    levels = c("normal", "heat"))
sample_data(ps3_order_bacteria)$condition <- factor(sample_data(ps3_order_bacteria)$condition, 
                                                    levels = c("normal", "heat"))
sample_data(ps4_order_bacteria)$condition <- factor(sample_data(ps4_order_bacteria)$condition, 
                                                    levels = c("normal", "heat", "extreme"))
sample_data(ps5_order_bacteria)$condition <- factor(sample_data(ps5_order_bacteria)$condition, 
                                                    levels = c("normal", "heat"))

## option 1 merge and run random-effects family-wise error correction 
ps_merged = merge_phyloseq(ps1_order_bacteria, ps2_order_bacteria, ps3_order_bacteria, ps4_order_bacteria, ps5_order_bacteria)
sample_data(ps_merged)$study
length(unique(sample_data(ps_merged)$study))

otu <- as.data.frame(otu_table(ps_merged))  # samples x taxa
otu$study <-  meta$study 
View(otu)


otu_psmerge <- as.data.frame(otu_table(ps_merged))
View(otu_psmerge)
study_coverage <- otu %>%
  group_by(study) %>%
  summarise(across(everything(), ~ sum(. > 0))) %>%
  column_to_rownames("study")
View(study_coverage) ## shows that aside from study 1, all studies have zero counts in the remaining ASVs
problematic_taxa <- names(which(colSums(study_coverage > 0) < 2))
length(problematic_taxa)

dim(study_coverage)

 ## option 2 run each independently and random effects meta 

sample_data(ps2_order_bacteria)
out <- ancombc2(
  data          = ps1_order_bacteria,
  assay_name    = "counts",
  tax_level     = "Family",        
  fix_formula   = "condition",
  rand_formula  = NULL,
  p_adj_method  = "BH",
  prv_cut       = 0.10,           
  lib_cut       = 1000,
  s0_perc       = 0.05,
  group         = "condition",
  struc_zero    = TRUE,
  neg_lb        = TRUE,
  alpha         = 0.05,
  n_cl          = 1
)
sample_data(ps2_order_bacteria)
out_2 <- ancombc2(
  data          = ps2_order_bacteria,
  assay_name    = "counts",
  tax_level     = "Family",        
  fix_formula   = "condition",
  rand_formula  = NULL,
  p_adj_method  = "BH",
  prv_cut       = 0.10,           
  lib_cut       = 1000,
  s0_perc       = 0.05,
  group         = "condition",
  struc_zero    = TRUE,
  neg_lb        = TRUE,
  alpha         = 0.05,
  n_cl          = 1
)

out_3 <- ancombc2(
  data          = ps3_order_bacteria,
  assay_name    = "counts",
  tax_level     = "Family",        
  fix_formula   = "condition + Species",
  rand_formula  = NULL,
  p_adj_method  = "BH",
  prv_cut       = 0.10,           
  lib_cut       = 1000,
  s0_perc       = 0.05,
  group         = "condition",
  struc_zero    = TRUE,
  neg_lb        = TRUE,
  alpha         = 0.05,
  n_cl          = 1
)


out_4 <- ancombc2(
  data          = ps4_order_bacteria,
  assay_name    = "counts",
  tax_level     = "Family",        
  fix_formula   = "condition + collection.time",
  rand_formula  = NULL,
  p_adj_method  = "BH",
  prv_cut       = 0.10,           
  lib_cut       = 1000,
  s0_perc       = 0.05,
  group         = "condition",
  struc_zero    = TRUE,
  neg_lb        = TRUE,
  alpha         = 0.05,
  n_cl          = 1
)

out_5 <- ancombc2(
  data          = ps5_order_bacteria,
  assay_name    = "counts",
  tax_level     = "Family",        
  fix_formula   = "condition + collection.time",
  rand_formula  = NULL,
  p_adj_method  = "BH",
  prv_cut       = 0.10,           
  lib_cut       = 1000,
  s0_perc       = 0.05,
  group         = "condition",
  struc_zero    = TRUE,
  neg_lb        = TRUE,
  alpha         = 0.05,
  n_cl          = 1
)

# save outputs for future ease 

saveRDS(out, 'out.rds')
saveRDS(out_2, 'out_2.rds')
saveRDS(out_3, 'out_3.rds')
saveRDS(out_4, 'out_4.rds')
saveRDS(out_5, 'out_5.rds')
out <- readRDS('out.rds')
out_2 <- readRDS('out_2.rds')
out_3 <- readRDS('out_3.rds')
out_4 <- readRDS('out_4.rds')
out_5 <- readRDS('out_5.rds')

# select for LFC columns that match the condition, in this case conditionheat
grep('^lfc_', names(out_3$res))
grep("^se_",  names(out$res))
grep("^q_",  names(out$res))
out$res

# function retrieves lfc, se and q values for each taxon associated with conditionheat
get_lfc <- function(res, study_name) {
  n <- names(res)
  data.frame(
    taxon = as.character(res$taxon),
    lfc   = as.numeric(res[[ n[grep("^lfc_conditionheat", n)]]]),
    se    = as.numeric(res[[ n[grep("^se_conditionheat", n)]]]),
    q     = as.numeric(res[[ n[grep("^q_conditionheat", n) ]]]),
    study = study_name
  )
}

# apply function to extract the lfc, se, and q estimates for the asvs 
lfc_ps1 <- get_lfc(out$res, "ps1")
lfc_ps2 <- get_lfc(out_2$res, "ps2")
lfc_ps3 <- get_lfc(out_3$res, "ps3")
lfc_ps4 <- get_lfc(out_4$res, "ps4")
lfc_ps5 <- get_lfc(out_5$res, "ps5")

View(lfc_ps5)
# check that the right condition was used 
res<- out$res
res$lfc_conditionnormal == lfc_ps1$lfc

# merge all of them together
lfc_long <- rbind(lfc_ps1, lfc_ps2, lfc_ps3, lfc_ps4, lfc_ps5)

dim(lfc_long)
head(lfc_long)
table(lfc_long$study)
view(lfc_long)

# select for ASVs seen in more than two groups 
shared_taxa <- lfc_long %>%
  filter(!is.na(lfc), se > 0) %>%
  count(taxon) %>%
  filter(n >= 3) %>%
  pull(taxon)

# filtering out any se <= 0 
lfc_filtered <- lfc_long %>%
  filter(taxon %in% shared_taxa, !is.na(lfc), q < 0.05)

lfc_long
lfc_filtered
shared_taxa
# Check which studies contribute to shared_taxa


data_lfc5 <- lfc_long %>% 
  filter(taxon %in% shared_taxa, !is.na(lfc), se > 0) %>%
  group_by(taxon) %>%
  summarise(studies = paste(study, collapse = ", "),
            n_studies = n()) %>% as.data.frame()

View(data_lfc5)

# random effects meta -analysis with moderator variables, using dplyr to keep unique taxons

taxa_order <- lfc_filtered %>%
  group_by(taxon) %>%
  group_keys() %>%
  pull(taxon)

meta_results <- lfc_filtered %>%
  group_by(taxon) %>%
  group_map(~ tryCatch(
    rma(yi = lfc, sei = se, data = .x, method = "REML"),
    error = function(e) NULL
  ), .keep = TRUE) %>%
  setNames(taxa_order)  


meta_summary <- imap_dfr(meta_results, function(m, taxon) {
  if (is.null(m)) return(NULL)
  tibble(
    taxon      = taxon,
    pooled_lfc = as.numeric(m$b),
    pooled_se  = as.numeric(m$se),
    ci_lo      = as.numeric(m$ci.lb),
    ci_hi      = as.numeric(m$ci.ub),
    pval       = as.numeric(m$pval),
    I2         = as.numeric(m$I2),
    k          = as.numeric(m$k)
  )
}) %>%
  mutate(q = p.adjust(pval, method = "BH")) %>%
  arrange(q, desc(abs(pooled_lfc))) %>%
  mutate(direction = case_when(
    q < 0.05 & pooled_lfc > 0 ~ "Enriched",
    q < 0.05 & pooled_lfc < 0 ~ "Depleted",
    TRUE                       ~ "Not significant"
  ))

View(meta_summary)

 
