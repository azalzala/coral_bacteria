library(dplyr)
library(ggplot2)
library(forcats)
library(phyloseq)
library(dplyr)
library(purrr)

ps_list <- list( "1" = ps1_order_bacteria, "2" = ps2_order_bacteria, "3" = ps3_order_bacteria, "4"=ps4_order_bacteria, "5" =ps5_order_bacteria)

asv_counts <- data.frame(
  dataset = names(ps_list),
  total_asvs = sapply(ps_list, ntaxa)
)
ggplot(asv_counts, aes(x = reorder(dataset, -total_asvs), y = total_asvs, fill = dataset)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = total_asvs), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Study",
    y = "Total ASVs",
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave("asv_counts.png", width = 8, height = 5, dpi = 300)


combo_summary <- lfc_long %>%
  filter(!is.na(lfc), se > 0) %>%
  group_by(taxon) %>%
  summarise(
    mean_lfc    = mean(lfc),
    n_studies   = n(),
    study_combo = paste(sort(study), collapse = " + ")
  ) %>%
  mutate(direction = ifelse(mean_lfc > 0, "Enriched", "Depleted"))

combo_counts <- combo_summary %>%
  group_by(study_combo, direction, n_studies) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    count_plot  = ifelse(direction == "Depleted", -count, count),
    combo_label = gsub("Study", "S", study_combo)  # shorten labels
  ) %>%
  arrange(n_studies, combo_label)

# Fix factor order for y-axis (grouped by n_studies)
combo_counts <- combo_counts %>%
  mutate(combo_label = fct_reorder(combo_label, n_studies))

# Divider positions between n_studies groups
dividers <- combo_counts %>%
  distinct(combo_label, n_studies) %>%
  arrange(n_studies, combo_label) %>%
  mutate(pos = as.numeric(as.factor(combo_label))) %>%
  group_by(n_studies) %>%
  summarise(max_pos = max(pos)) %>%
  filter(n_studies < max(n_studies)) %>%
  pull(max_pos)

# Group midpoint annotations
group_mids <- combo_counts %>%
  distinct(combo_label, n_studies) %>%
  arrange(n_studies) %>%
  mutate(pos = as.numeric(fct_reorder(combo_label, n_studies))) %>%
  group_by(n_studies) %>%
  summarise(mid = mean(pos))

ggplot(combo_counts, aes(x = count_plot, y = combo_label, fill = direction)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.8, colour = "black") +
  
  # Dividers between study-count groups
  geom_hline(yintercept = dividers + 0.5,
             linetype = "dotted", colour = "grey60", linewidth = 0.7) +
  
  # n_studies annotations on right
  geom_text(data = group_mids,
            aes(x = Inf, y = mid, label = paste0(n_studies, " studies")),
            inherit.aes = FALSE,
            hjust = -0.1, size = 3.2, colour = "grey40", fontface = "bold") +
  
  scale_fill_manual(values = c("Enriched" = "#4169E1", "Depleted" = "#CD7054")) +
  scale_x_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.05, 0.2))  # room for annotations
  ) +
  labs(
    x        = "Number of Taxa",
    y        = "Study Combination",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(colour = "grey40", size = 10),
    axis.text.y      = element_text(size = 10),
    plot.margin      = margin(10, 60, 10, 10)  # right margin for annotations
  )

ggsave("study_combination_taxons.png", width = 10, height = 7, dpi = 300)

sum(unique(lfc_long$taxon) %in% tax_df$taxon)  # should equal 
length(unique(lfc_long$taxon))
length(unique(tax_df$taxon))

tax_map <- tax_df %>%
  select(taxon, Family Genus) %>%
  deframe()

write.csv(tax_map, 'tax_map.csv')

View(tax_df)

# In R, export your taxonomy
tax_df <- as.data.frame(tax_table(ps_merged)) %>%
  rownames_to_column("taxon") %>%
  select(taxon, Family, Genus)

write.csv(tax_df, "tax_map.csv", row.names = FALSE)
write.csv(lfc_long, "lfc_long.csv", row.names = FALSE)
head(tax_df$taxon)
head(lfc_long$taxon)

tax_df <- as.data.frame(tax_table(ps_merged)) %>%
  rownames_to_column("taxon") %>%   # rownames = whatever lfc_long$taxon uses
  select(taxon, Family, Genus)

write.csv(tax_df, "tax_map.csv", row.names = FALSE)
sum(unique(lfc_long$taxon) %in% tax_df$taxon)

tax_df <- as.data.frame(tax_table(ps_merged)) %>%
  rownames_to_column("sequence") %>%
  mutate(taxon = Genus) %>%
  select(taxon, Family, Genus) %>%
  distinct(taxon, .keep_all = TRUE)# deduplicate since multiple ASVs per genus

write.csv(tax_df, "tax_map.csv", row.names = FALSE)

# shared_taxa needs these columns
shared_taxa_export <- meta_summary %>%
  filter(taxon %in% shared_taxa) %>%
  select(taxon, pooled_lfc, ci_lo, ci_hi, q, k) %>%
  left_join(
    lfc_long %>%
      filter(!is.na(lfc), se > 0) %>%
      group_by(taxon) %>%
      summarise(studies   = paste(sort(study), collapse = " + "),
                n_studies = n()),
    by = "taxon"
  )

write.csv(shared_taxa_export, "shared_taxa.csv", row.names = FALSE)

view(shared_taxa_export)


# SIGNIFICANT TAXON ENRICHED/DEPLETED FOREST PLOT
ggplot(meta_summary, aes(x = pooled_lfc, y = -log10(q), colour = direction)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "black") +
  scale_colour_manual(
    values = c("Enriched" = "#CD7054", "Depleted" = "#4169E1", "Not significant" = "grey70"),
    name = "Direction"
  ) +
  labs(
    x = "Pooled Log-Fold Change",
    y = expression(-log[10](q)),
    title = "Meta-Analysis Differential Abundance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("waterfall_plot.png", width = 11, height = 6, dpi = 300)


waterfall_data <- meta_summary %>%
  filter(taxon %in% shared_taxa, q < 0.05) %>%
  arrange(pooled_lfc) %>%
  mutate(direction = ifelse(pooled_lfc > 0, "Enriched", "Depleted"))

# Lock factor order separately after arrange
waterfall_data$taxon <- factor(waterfall_data$taxon, 
                               levels = waterfall_data$taxon)

nrow(waterfall_data)          # should be > 0
levels(waterfall_data$taxon)  # should show taxa in LFC order
table(waterfall_data$direction)

ggplot(waterfall_data, aes(x = taxon, y = pooled_lfc, fill = direction)) +
  
  # Bars
  geom_col(width = 0.7) +
  
  # 95% CI error bars
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.25, linewidth = 0.5, colour = "grey30") +
  
  # Zero line
  geom_hline(yintercept = 0, linewidth = 0.8, colour = "black") +
  
  # Significance stars
  geom_text(aes(
    label = case_when(q < 0.001 ~ "***",
                      q < 0.01  ~ "**",
                      q < 0.05  ~ "*",
                      TRUE      ~ ""),
    y = ifelse(pooled_lfc > 0, ci_hi + 0.1, ci_lo - 0.1)
  ), size = 4, vjust = 0.5) +
  
  scale_fill_manual(values = c("Enriched" = "#4169E1", "Depleted" = "#CD7054")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  
  labs(
    x        = NULL,
    y        = "Pooled Log Fold Change (LFC)",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    legend.position  = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(colour = "grey40", size = 10)
  )

# ROBUST MARKERS BAR PLOT

meta_summary_qfilter <- meta_summary %>%
  filter(q < 0.05, I2 < 50)

robust_markers <- ggplot(meta_summary_qfilter, 
                         aes(x = reorder(taxon, pooled_lfc), 
                             y = pooled_lfc, 
                             fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Enriched" = "#CD7054", "Depleted" = "#4169E1"),
    name = "Direction"
  ) +
  labs(
    x = "Species",
    y = "Pooled LFC Estimate"
  ) +
  theme_classic()

robust_markers

ggsave('robust_markers.png', plot = robust_markers )

## PACKAGES FOR HEATMAP
library(ggplot2)
library(dplyr)
library(tidyr)

# filtering for statistically significant taxa and selecting taxon
sig_taxa <- meta_summary %>%
  filter(q < 0.05) %>%
  pull(taxon)

# Pivot to wide format: taxa × studies
heatmap_data <- lfc_filtered %>%
  filter(taxon %in% sig_taxa) %>%
  select(taxon, study, lfc) %>%
  pivot_wider(names_from = study, values_from = lfc)

# Back to long for ggplot
heatmap_long <- heatmap_data %>%
  pivot_longer(-taxon, names_to = "study", values_to = "lfc")

# Clamp extreme LFC values for colour scale
heatmap_long <- heatmap_long %>%
  mutate(lfc_clamped = pmax(pmin(lfc, 3), -3))

ggplot(heatmap_long, aes(x = study, y = reorder(taxon, lfc, mean, na.rm = TRUE), 
                         fill = lfc_clamped)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low  = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey90",
    name = "LFC\n(clamped ±3)"
  ) +
  labs(
    title    = "Log-fold change of coral microbiome taxa under heat stress",
    subtitle = "Significant taxa (q < 0.05) across studies",
    x        = "Study",
    y        = "Taxon"
  ) +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 7),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank()
  )

## BAR PLOT TAXONOMY

# Aggregate to a chosen taxonomic level (e.g. Phylum or Family)
tax_level <- "Phylum"  # change to "Family", "Genus" etc.

# Function to extract relative abundance per study
# Aggregate to relative abundance per study, grouped by condition
get_rel_abund <- function(ps, study_name) {
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  
  psmelt(ps_rel) %>%
    group_by(condition, .data[[tax_level]]) %>%  # ← group by condition, not Sample
    summarise(Abundance = mean(Abundance), .groups = "drop") %>%
    mutate(study = study_name)
}

# Apply across all studies
ps_list <- list("study 1" = ps1, "Study 2" = ps2, "study 3" = ps3,
                "study 4" = ps4, "study 5" = ps5)

bar_data <- imap_dfr(ps_list, ~ get_rel_abund(.x, .y))

top_taxa <- 12  # ← updated

bar_data <- bar_data %>%
  mutate(!!tax_level := ifelse(.data[[tax_level]] %in% top_taxa,
                               .data[[tax_level]], "Other"))

ggplot(bar_data, aes(x = condition, y = Abundance, fill = .data[[tax_level]])) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ study, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Paired", na.value = "grey80") +
  labs(
    title = paste("Taxonomic composition by study —", tax_level, "level"),
    x     = "Condition",
    y     = "Relative Abundance",
    fill  = tax_level
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold"),
    legend.position = "bottom"
  )

