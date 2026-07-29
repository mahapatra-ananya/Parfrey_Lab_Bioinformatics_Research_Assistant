library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ALDEx2)
library(dplyr)
library(ggdendro)
library(ggrepel)
library(phyloseq)
library(MicrobiomeStat)
# library(Maaslin2)
library(lefser)
library(forcats)
library(qualpalr)
library("DESeq2")
library(vegan)

# First read in and subset abundance dataframes and metadata to only rhizomes on day 20 with the 4 treatments we care about
subset_naming_pattern <- c("function", "pathway", "description", "_SI", "_SR", "_CI", "_CR")
subset_naming_pattern_sulfidic <- c("function", "pathway", "description", "_SI", "_SR")
subset_naming_pattern_control <- c("function", "pathway", "description", "_CI", "_CR")

# full_metacyc <- read_delim("Results/picrust2_output/pathways_out/path_abun_unstrat_descrip.tsv", delim="\t", col_names=TRUE, trim_ws=TRUE)
# subset_metacyc <- full_metacyc %>% 
#   dplyr::select(ends_with(subset_naming_pattern))
# 
# full_ko <- read_delim("Results/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv", delim="\t", col_names=TRUE, trim_ws=TRUE)
# subset_ko <- full_ko %>% 
#   dplyr::select(ends_with(subset_naming_pattern))

full_ec <- read_delim("Results/picrust2_output/ec_metagenome_out/pred_metagenome_unstrat.tsv", delim="\t", col_names=TRUE, trim_ws=TRUE)
subset_ec <- full_ec %>% 
  dplyr::select(ends_with(subset_naming_pattern))
subset_ec_sulfidic <- full_ec %>% 
  dplyr::select(ends_with(subset_naming_pattern_sulfidic))
subset_ec_control <- full_ec %>% 
  dplyr::select(ends_with(subset_naming_pattern_control))

# full_taxa_contrib <- read_delim("Results/picrust2_output/ec_metagenome_out/pred_metagenome_contrib.tsv", delim="\t", col_names=TRUE, trim_ws=TRUE)
# suffix_pattern <- paste0("(", paste(c("_SI", "_SR", "_CI", "_CR"), collapse = "|"), ")$")
# subset_taxa_contrib <- full_taxa_contrib %>% 
#   dplyr::filter(grepl(suffix_pattern, sample))

# write_tsv(subset_metacyc, "Results/picrust2_output/pathways_out/path_abun_unstrat_descrip_subset.tsv")
# write_tsv(subset_ko, "Results/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat_descrip_subset.tsv")
# write_tsv(subset_ko, "Results/picrust2_output/EC_metagenome_out/pred_metagenome_unstrat_descrip_subset.tsv")
# 
# metacyc_file <- "Results/picrust2_output/pathways_out/path_abun_unstrat_descrip_subset.tsv" # file path
# ko_file <- "Results/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat_descrip_subset.tsv" # file path
# ec_file <- "Results/picrust2_output/EC_metagenome_out/pred_metagenome_unstrat_descrip_subset.tsv" # file path

metadata <- read_delim(
  "metadata_phyloseq.tsv",
  delim = "\t",
  escape_double = FALSE, # doesn't matter for our metadata
  trim_ws = TRUE
)
subset_metadata <- metadata %>% filter(sample_type == "rhizome" & experiment_day == 20)
subset_metadata_sulfidic <- subset_metadata %>% filter(sed_treat == "sulfidic")
subset_metadata_control <- subset_metadata %>% filter(sed_treat == "control")

subset_ec <- pathway_annotation(
  pathway = "EC",
  daa_results_df = subset_ec,
  ko_to_kegg = FALSE
)
subset_ec_sulfidic <- pathway_annotation(
  pathway = "EC",
  daa_results_df = subset_ec_sulfidic,
  ko_to_kegg = FALSE
)
subset_ec_control <- pathway_annotation(
  pathway = "EC",
  daa_results_df = subset_ec_control,
  ko_to_kegg = FALSE
)

##################### Looking for Julia's EC pathways in picrust2 output ###############
## 2026-07-27 update: includes iron pathways 

full_relevant_ec_pathways_df <- read_delim(
  "Results/picrust2_output/from_julia/matched_genes.tsv",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
) 

iron_pathways_df <- read_delim(
  "Results/picrust2_output/from_julia/iron_pathways_genes.tsv",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_names = FALSE
)
iron_pathways_df <- iron_pathways_df[,c(1,3)]
colnames(iron_pathways_df) <- c("pathway_name", "ec_number")

relevant_ec_pathways_df <- full_relevant_ec_pathways_df[, c(1,3)] # just keep pathway name and ec number
# many rows will be the same, as the full df contains different genes for the same pathway and ec number
# Total number of duplicate rows:
relevant_ec_pathways_df <- rbind(relevant_ec_pathways_df, iron_pathways_df)
sum(duplicated(relevant_ec_pathways_df))
relevant_ec_pathways_df <- unique(relevant_ec_pathways_df) # remove exact duplicate rows
# now some rows will have the same ec number because it maps to different pathways
# change pathway_name to be each pathway separated by commas in the same cell
relevant_ec_pathways_df <- relevant_ec_pathways_df %>%
  group_by(ec_number) %>%
  summarise(pathway_name = paste0(pathway_name, collapse = ", "), .groups = "drop")

relevant_ec_picrust_df <- subset_ec %>% 
  dplyr::rename("feature" = "function") %>% 
  filter(feature %in% relevant_ec_pathways_df[["ec_number"]])

########################## Running ggpicrust2 analysis ########################## 

run_ggpicrust2 <- function(df, metadata, group, method, reference) {
  results <- ggpicrust2(data = df,
                        metadata = metadata,
                        group = group,
                        pathway = "EC",
                        daa_method = method,
                        ko_to_kegg = FALSE,
                        order = "group",
                        reference = reference,
                        p_values_bar = TRUE,
                        x_lab = "description",
                        filter_for_prokaryotes = TRUE)
  return(results)
}

relevant_ec_results_aldex2 <- run_ggpicrust2(relevant_ec_picrust_df, subset_metadata, "sed_treat", "ALDEx2", "control") 
relevant_ec_daa_results_df_welch <- subset(relevant_ec_results_aldex2[["daa_results_df"]], method == "ALDEx2_Welch's t test")
relevant_ec_daa_results_df_wilcox <- subset(relevant_ec_results_aldex2[["daa_results_df"]], method == "ALDEx2_Wilcoxon rank test")
relevant_ec_results_deseq2 <- run_ggpicrust2(relevant_ec_picrust_df, subset_metadata, "sed_treat", "DESeq2", "control")
relevant_ec_results_linda <- run_ggpicrust2(relevant_ec_picrust_df, subset_metadata, "sed_treat", "LinDA", "control")

relevant_ec_pathways_df <- relevant_ec_pathways_df %>%
  dplyr::rename("feature" = "ec_number")
linda_daa_results_df <- relevant_ec_results_linda[["daa_results_df"]] %>%
  inner_join(relevant_ec_pathways_df)
aldex2_welch_daa_results_df <- relevant_ec_daa_results_df_welch %>%
  inner_join(relevant_ec_pathways_df)
aldex2_wilcox_daa_results_df <- relevant_ec_daa_results_df_wilcox %>%
  inner_join(relevant_ec_pathways_df)
deseq2_daa_results_df <- relevant_ec_results_deseq2[["daa_results_df"]] %>%
  inner_join(relevant_ec_pathways_df)

combined_daa_results_df <- linda_daa_results_df %>%
  inner_join(deseq2_daa_results_df, by = "feature") %>%
  inner_join(aldex2_welch_daa_results_df, by = "feature") %>%
  inner_join(aldex2_wilcox_daa_results_df, by = "feature") %>%
  dplyr::select(
    "ec_number" = "feature",
    "description" = "description.x",
    "pathway_name" = "pathway_name.x",
    "group1" = "group1.x",
    "group2" = "group2.x",
    "l2fc_LinDA" = "log2_fold_change.x",
    "pval_LinDA" = "p_values.x",
    "padj_LinDA" = "p_adjust.x",
    "l2fc_DESeq2" = "log2_fold_change.y",
    "pval_DESeq2" = "p_values.y",
    "padj_DESeq2" = "p_adjust.y",
    "l2fc_ALDEx2_Welch" = "log2_fold_change.x.x",
    "pval_ALDEx2_Welch" = "p_values.x.x",
    "padj_ALDEx2_Welch" = "p_adjust.x.x",
    "l2fc_ALDEx2_Wilcox" = "log2_fold_change.y.y",
    "pval_ALDEx2_Wilcox" = "p_values.y.y",
    "padj_ALDEx2_Wilcox" = "p_adjust.y.y"
  )

write.csv(linda_daa_results_df, "linda_daa_results.csv")
write.csv(deseq2_daa_results_df, "deseq2_daa_results.csv")
write.csv(aldex2_welch_daa_results_df, "aldex2_welch_daa_results.csv")
write.csv(aldex2_wilcox_daa_results_df, "aldex2_wilcox_daa_results.csv")
write.csv(combined_daa_results_df, "combined_daa_results.csv")

########################## Comparing DA methods ########################## 
methods <- c("LinDA", "DESeq2", "ALDEx2")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = relevant_ec_picrust_df %>% column_to_rownames("description"), 
              metadata = subset_metadata, group = "sed_treat", 
              daa_method = method, reference = "control")
})
method_names <- c("LinDA", "DESeq2", "ALDEx2")
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = method_names)

aldex2_welch_sig_df <- subset(daa_results_list[[3]], p_adjust <= 0.05 & method == "ALDEx2_Welch's t test")
aldex2_wilcox_sig_df <- subset(daa_results_list[[3]], p_adjust <= 0.05 & method == "ALDEx2_Wilcoxon rank test")
deseq2_sig_df <- subset(daa_results_list[[2]], p_adjust <= 0.05)
linda_sig_df <- subset(daa_results_list[[1]], p_adjust <= 0.05)

common_sig_df <- linda_sig_df %>%
  inner_join(deseq2_sig_df, by="feature") %>%
  inner_join(aldex2_welch_sig_df, by="feature") %>%
  inner_join(aldex2_wilcox_sig_df, by="feature") %>%
  dplyr::select(
    "description" = "feature",
    "group1" = "group1.x",
    "group2" = "group2.x",
    "l2fc_LinDA" = "log2_fold_change.x",
    "pval_LinDA" = "p_values.x",
    "padj_LinDA" = "p_adjust.x",
    "l2fc_DESeq2" = "log2_fold_change.y",
    "pval_DESeq2" = "p_values.y",
    "padj_DESeq2" = "p_adjust.y",
    "l2fc_ALDEx2_Welch" = "log2_fold_change.x.x",
    "pval_ALDEx2_Welch" = "p_values.x.x",
    "padj_ALDEx2_Welch" = "p_adjust.x.x",
    "l2fc_ALDEx2_Wilcox" = "log2_fold_change.y.y",
    "pval_ALDEx2_Wilcox" = "p_values.y.y",
    "padj_ALDEx2_Wilcox" = "p_adjust.y.y"
  ) %>%
  mutate(group1 = "control", group2 = "sulfidic") %>%
  inner_join(combined_daa_results_df %>% dplyr::select("ec_number", "description", "pathway_name"), by = "description") %>%
  relocate("ec_number", "pathway_name")

write.csv(common_sig_df, "common_sig_pathways_df.csv")

combined_sig_df <- linda_sig_df %>%
  full_join(deseq2_sig_df, by="feature") %>%
  full_join(aldex2_welch_sig_df, by="feature") %>%
  full_join(aldex2_wilcox_sig_df, by="feature") %>%
  dplyr::select(
    "description" = "feature",
    "group1" = "group1.x",
    "group2" = "group2.x",
    "l2fc_LinDA" = "log2_fold_change.x",
    "pval_LinDA" = "p_values.x",
    "padj_LinDA" = "p_adjust.x",
    "l2fc_DESeq2" = "log2_fold_change.y",
    "pval_DESeq2" = "p_values.y",
    "padj_DESeq2" = "p_adjust.y",
    "l2fc_ALDEx2_Welch" = "log2_fold_change.x.x",
    "pval_ALDEx2_Welch" = "p_values.x.x",
    "padj_ALDEx2_Welch" = "p_adjust.x.x",
    "l2fc_ALDEx2_Wilcox" = "log2_fold_change.y.y",
    "pval_ALDEx2_Wilcox" = "p_values.y.y",
    "padj_ALDEx2_Wilcox" = "p_adjust.y.y"
  ) %>%
  mutate(group1 = "control", group2 = "sulfidic") %>%
  inner_join(combined_daa_results_df %>% dplyr::select("ec_number", "description", "pathway_name"), by = "description") %>%
  relocate("ec_number", "pathway_name")

write.csv(combined_sig_df, "combined_sig_pathways_df.csv")

########################## Visualizing DA ########################## 

## PCA
relevant_ec_pca_plot <- pathway_pca(
  abundance = relevant_ec_results_aldex2$abundance,
  metadata = relevant_ec_results_aldex2$metadata,
  group = relevant_ec_results_aldex2$group
)
relevant_ec_pca_plot

## Functions to visualize DA per method

make_heatmap <- function(sig_df, sig_abund, metadata) {
  if (length(sig_df$feature) > 0) {
    result <- pathway_heatmap(
      abundance = sig_abund,
      metadata = metadata,
      group = "treat",
      cluster_rows = FALSE, # cluster according to our arrangement!! not internally by abund
      cluster_cols = TRUE # dendrogram
    )
  }
  label_mapping <- setNames(sig_df$combined_label, sig_df$feature)
  result <- result + scale_y_discrete(labels = label_mapping)
  return(result)
}

make_volcanoplot <- function(df, da_method) {
  result <- pathway_volcano(
    df,
    fc_col = "log2_fold_change",
    p_col = "p_adjust",
    label_col = "combined_label",
    p_threshold = 0.05,
    fc_threshold = 0.5, # can change to 1, using 0.5 to see what else pops up
    label_top_n = 40,
    point_size = 2,
    point_alpha = 0.6,
    colors = c(Down = "blue", `Not Significant` = "grey60", Up = "red"),
    show_threshold_lines = TRUE,
    title = "Volcano Plot: Relevant EC Pathway Differential Abundance",
    x_lab = "log2 Fold Change",
    y_lab = paste0("-log10(Adjusted P-value) ", da_method)
  )
  return(result)
}

make_errorbarplot <- function(df, abundance, group) {
  error_df <- df %>% mutate(description = combined_label)
  plot <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = error_df,        
    Group = group,         
    p_values_threshold = 0.05,
    order = "p_values",                         
    p_value_bar = TRUE,
    colors = NULL,                              
    x_lab = "description"               
  )
  return(plot)
}

make_plots <- function(df, metadata, da_method, ggpicrust_results, group) {
  df <- df %>%
    left_join(
      relevant_ec_pathways_df,            
      by = "feature" 
    ) %>%
    dplyr::arrange(pathway_name, p_adjust) %>% # arranging by pathway name
    mutate(combined_label = paste0(pathway_name, ": ", description)) # combined label
  sig_df <- df %>%
    dplyr::filter(p_adjust < 0.05)  # filter first 
  # dplyr::arrange(p_adjust)        # rank within significant set
  # dplyr::slice_head(n = 30)           # take top 30 most significant
  sig_abund <- ggpicrust_results$abundance[sig_df$feature, , drop = FALSE]
  heatmap <- make_heatmap(sig_df, sig_abund, ggpicrust_results$metadata)
  volcanoplot <- make_volcanoplot(df, da_method)
  errorbarplot <- make_errorbarplot(df, ggpicrust_results$abundance, metadata[[group]])
  return(list("heatmap" = heatmap, "volcanoplot" = volcanoplot, "errorbarplot" = errorbarplot))
}

plots_welch <- make_plots(relevant_ec_daa_results_df_welch, subset_metadata, "ALDEx2 Welch", relevant_ec_results_aldex2, "sed_treat")
plots_welch$heatmap
plots_welch$volcanoplot
plots_welch$errorbarplot

plots_wilcox <- make_plots(relevant_ec_daa_results_df_wilcox, subset_metadata, "ALDEx2 Wilcox", relevant_ec_results_aldex2, "sed_treat")
plots_wilcox$heatmap
plots_wilcox$volcanoplot
plots_wilcox$errorbarplot

plots_deseq2 <- make_plots(relevant_ec_results_deseq2[["daa_results_df"]], subset_metadata, "DESeq2", relevant_ec_results_deseq2, "sed_treat")
plots_deseq2$heatmap
plots_deseq2$volcanoplot
plots_deseq2$errorbarplot

plots_linda <- make_plots(relevant_ec_results_linda[["daa_results_df"]], subset_metadata, "LinDA", relevant_ec_results_linda, "sed_treat")
plots_linda$heatmap
plots_linda$volcanoplot
plots_linda$errorbarplot



# saveRDS(relevant_ec_results, "relevant_ec_ggpicrust_results.RDS")
# saveRDS(ec_results, "ec_ggpicrust_results.RDS")
# saveRDS(ko_results, "ko_ggpicrust_results.RDS")

######################################## sulfidic ############################ 

relevant_ec_picrust_df_sulf <- subset_ec_sulfidic %>% 
  dplyr::rename("feature" = "function") %>%
  filter(feature %in% relevant_ec_pathways_df[["feature"]])

relevant_ec_results_sulfidic <- run_ggpicrust2(relevant_ec_picrust_df_sulf, 
                                               subset_metadata_sulfidic, "microbe_treat", 
                                               "LinDA", "intact")

## PCA
relevant_ec_pca_plot_sulf <- pathway_pca(
  abundance = relevant_ec_results_sulfidic$abundance,
  metadata = relevant_ec_results_sulfidic$metadata,
  group = relevant_ec_results_sulfidic$group
)
relevant_ec_pca_plot_sulf

plots_sulfidic <- make_plots(relevant_ec_results_sulfidic[["daa_results_df"]], subset_metadata_sulfidic, "LinDA", relevant_ec_results_sulfidic, "microbe_treat")
plots_sulfidic$heatmap
plots_sulfidic$volcanoplot
plots_sulfidic$errorbarplot

######################################## control ############################ HAVE TO MAKE FUNCTION ASAP

relevant_ec_picrust_df_control <- subset_ec_control %>% 
  dplyr::rename("feature" = "function") %>% 
  filter(feature %in% relevant_ec_pathways_df[["feature"]])

relevant_ec_results_control <- run_ggpicrust2(relevant_ec_picrust_df_control, 
                                               subset_metadata_control, "microbe_treat", 
                                               "DESeq2", "intact")

## PCA
relevant_ec_pca_plot_control <- pathway_pca(
  abundance = relevant_ec_results_control$abundance,
  metadata = relevant_ec_results_control$metadata,
  group = relevant_ec_results_control$group
)
relevant_ec_pca_plot_control

plots_control <- make_plots(relevant_ec_results_control[["daa_results_df"]], subset_metadata_control, "DESeq2", relevant_ec_results_control, "microbe_treat")
plots_control$heatmap
plots_control$volcanoplot
plots_control$errorbarplot



######################################## TAXA CONTRIB ###########################################

ps <- readRDS("../../Data/filtered_seagrass_mesocosm.rds")
ps <- subset_samples(ps, sample_type == "rhizome" & experiment_day == 20)
taxa_df <- as.data.frame(tax_table(ps))
taxa_df$taxon <- rownames(taxa_df)

# From picrust2 tutorial:
# taxon_abun - Abundance of this taxon in the sample. If abundances were normalized 
# by marker gene abundance this will be the normalized abundance, but NOT in terms of relative abundance.
# 
# taxon_rel_abun - This is the same as the "taxon_abun" column, but in terms of 
# relative abundance (so that the sum of all taxa abundances per sample is 100).
# 
# genome_function_count - Predicted copy number of this function per taxon.
# 
# taxon_function_abun - Multiplication of "taxon_abun" column by "genome_function_count" column.
# 
# taxon_rel_function_abun - Multiplication of "taxon_rel_abun" column by "genome_function_count" column.
# 
# norm_taxon_function_contrib - Proportional relative abundance of the taxon_function_abun 
# column per function and sample. I.e., it is what proportion of the specified function is 
# contributed by that particular taxon in the sample.


full_taxa_contrib <- read_contrib_file(file = "Results/picrust2_output/ec_metagenome_out/pred_metagenome_contrib.tsv")
suffix_pattern <- paste0("(", paste(c("_SI", "_SR", "_CI", "_CR"), collapse = "|"), ")$")
subset_taxa_contrib <- full_taxa_contrib %>% 
  dplyr::filter(grepl(suffix_pattern, sample))
aggregated_taxa_df <- aggregate_taxa_contributions(
  subset_taxa_contrib,
  #taxonomy = taxa_df,
  #tax_level = "Family",
  top_n = 100000000, # 107 to get 50 families
  daa_results_df = relevant_ec_results_linda[["daa_results_df"]],
  p_threshold = 0.05,
  contribution_col = "norm_taxon_function_contrib"
)

relevant_ec_picrust_df <- relevant_ec_picrust_df %>%
  left_join(
    relevant_ec_pathways_df,            
    by = "feature" 
  ) %>%
  dplyr::arrange(pathway_name) %>%
  mutate(combined_label = paste0(pathway_name, ": ", description))

family_map <- taxa_df %>%
  rownames_to_column(var = "taxon_id") %>%
  dplyr::select(taxon_id, Family)

# Match and replace the labels in aggregated_taxa_df
aggregated_taxa_df <- aggregated_taxa_df %>%
  # Left join to bring the family column in 
  left_join(family_map, by = c("taxon_label" = "taxon_id")) %>%
  # If a family was found, use it otherwise keep the original label ("Other" / "unclassified")
  mutate(taxon_label = ifelse(!is.na(Family), Family, taxon_label)) %>%
  # Clean up the extra column
  dplyr::select(-Family)

pathway_map <- relevant_ec_picrust_df %>%
  dplyr::select(feature, pathway_name, description, combined_label)

# Join and add the new column to aggregated dataframe
aggregated_taxa_df_mapped <- aggregated_taxa_df %>%
  left_join(pathway_map, by = c("function_id" = "feature")) %>%
  mutate(function_id = description) %>%
  dplyr::arrange(pathway_name)

aggregated_taxa_df_to_save <- aggregated_taxa_df %>%
  dplyr::select(
    "ec_number" = "function_id",
    "pathway_name",
    "description",
    "sample",
    "taxon family" = "taxon",
    "proportional_contrib_in_sample" = "contribution"
  )

write.csv(aggregated_taxa_df_to_save, "proportional_taxon_contrib.csv")

taxa_contribution_bar(
  contrib_agg = aggregated_taxa_df_mapped,
  metadata = subset_metadata,
  group = "sed_treat",
  facet_by = "function",
  n_functions = 50
)

taxa_contribution_heatmap(
  contrib_agg = aggregated_taxa_df_mapped,
  annotation_data = NULL,
  n_functions = 50,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  # clustering_method = "complete",
  # clustering_distance = "euclidean",
  low_color = "lightyellow",
  high_color = "green4",
  custom_title = "Proportional Taxa Contributions to Significant Relevant EC Pathways"
)

## Subset to sulfur cycling pathways only

aggregated_taxa_df <- aggregated_taxa_df %>%
  inner_join(pathway_map, by = c("function_id" = "feature"))

sc_aggregated_taxa_df <- aggregated_taxa_df %>%
  dplyr::filter(grepl("sulf", pathway_name, ignore.case = TRUE))

top_taxa <- sc_aggregated_taxa_df %>%
  group_by(combined_label, taxon) %>%
  summarise(
    avg_contribution = mean(contribution, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(combined_label) %>%
  slice_max(avg_contribution, n = 5, with_ties = TRUE) %>%
  dplyr::select(combined_label, taxon)

highest_contribution_df <- sc_aggregated_taxa_df %>%
  semi_join(
    top_taxa,
    by = c("combined_label", "taxon")
  )

plot_highest_contrib_df <- highest_contribution_df %>%
  group_by(combined_label, taxon) %>%
  summarise(
    avg_contribution = mean(contribution, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- plot_highest_contrib_df %>%
  mutate(
    taxon = as.factor(taxon)
  ) %>%
  group_by(combined_label) %>%
  mutate(
    taxon = fct_reorder(taxon, avg_contribution, .desc = TRUE)
  ) %>%
  ungroup() %>%
  mutate(taxon = fct_drop(taxon))

ggplot(plot_df, aes(
  x = taxon,
  y = avg_contribution
)) +
  geom_col() +
  facet_wrap(~ combined_label, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "Taxon", y = "Avg contribution")

ps_family <- ps
tax_df <- as.data.frame(tax_table(ps_family))

# define families to keep (these are family names)
selected_taxa <- unique(plot_highest_contrib_df$taxon)

# create grouping variable at family level
tax_df$Family <- ifelse(
  tax_df$Family %in% selected_taxa,
  tax_df$Family,
  "Others"
)

tax_table(ps_family) <- as.matrix(tax_df)

ps_grouped <- tax_glom(ps_family, taxrank = "Family", NArm = FALSE)

ps_rel <- transform_sample_counts(ps_grouped, function(x) x / sum(x))
plot_abund_df <- psmelt(ps_rel)
plot_abund_df$Family <- factor(
  plot_abund_df$Family,
  levels = c(setdiff(unique(plot_abund_df$Family), "Others"), "Others")
)

numcol <- length(unique(plot_abund_df$Family))
set.seed(15)
hex = as.data.frame(c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
           "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
           "#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363",
           "#a1d99b", "#fdae6b", "#9ecae1", "#fdd0a2", "#c7e9c0",
           "#f7f7f7", "#5254a3", "#6b6ecf", "#9c9ede", "#d6616b",
           "#ce6dbd", "#de9ed6", "#31a354"
))
colnames(hex) <- c("taxa_color")
tops = as.data.frame(c(unique(plot_abund_df$Family)))
colnames(tops) <- c("plotnames")
topcolors = cbind(tops, hex)
topcolors[topcolors$plotnames == "Others",]$taxa_color <- "grey90"
plotcolors <- topcolors$taxa_color
names(plotcolors) <- topcolors$plotnames

ggplot(plot_abund_df, aes(
  x = sample_name,
  y = as.numeric(Abundance),
  fill = Family
)) +
  geom_bar(stat = "identity") +
  
  # ONLY force grey for Others; everything else auto-colors
  scale_fill_manual(
    values = plotcolors
  ) +
  # 
  guides(fill = guide_legend(ncol = 2)) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(color = "black", size = 10),
    legend.text = element_text(size = 6),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank()
  ) +
  labs(
    y = "Relative Abundance",
    x = "Sample",
    fill = "Taxa"
  ) +
  facet_grid(. ~ treat, scales = "free")


# deseq2 DA on taxa contrib to sulfur cycling pathways
dds = phyloseq_to_deseq2(ps_grouped, ~ sed_treat)
# Run the actual differential abundance analysis using Wald test and positive counts (useful for sparse data)
dds = DESeq(dds, test="Wald", fitType="parametric", sfType="poscounts")

# Extract results
# cooksCutoff FALSE so it doesn't filter outliers automatically
# contrast explicitly defines numerator and denominator for easy interpretation of logfold2 change
res = results(dds, cooksCutoff = FALSE, contrast=c("sed_treat", "control", "sulfidic"))
# Table of significant results 
sigtab = res[which(res$padj < 0.05), ]
# Merge with taxonomy table (cbind appends columns)
res = cbind(as(res, "data.frame"), as(tax_table(ps_grouped)[rownames(res), ], "matrix"))
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_grouped)[rownames(sigtab), ], "matrix"))

# Create df for volcano plot
res <- as.data.frame(res) # convert to df
volcano_df <- res %>%
  filter(!is.na(padj))

volcano_df$Significant <- "Not Significant"
volcano_df$Significant[volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 0] <- "Control Enriched"
volcano_df$Significant[volcano_df$padj < 0.05 & volcano_df$log2FoldChange < 0] <- "Sulfidic Enriched"

# Below is copied from MJ's code for volcano plot from DESeq2 results
# In her 4_microbialanalysis.R file

# Extract OTU table
rel_mat <- as.data.frame(as(otu_table(ps_rel), "matrix"))

# Ensure taxa are rows
if (!taxa_are_rows(ps_rel)) {
  rel_mat <- t(rel_mat)
}

# Calculate mean relative abundance per ASV
mean_rel_abund <- rowMeans(rel_mat)

# Convert to dataframe for merging
rel_df <- data.frame(
  rel_abundance = mean_rel_abund
)

# Merge into volcano_df
#volcano_df <- merge(volcano_df, rel_df, by = "ASV", all.x = TRUE)
volcano_df = cbind(volcano_df, rel_df[rownames(volcano_df), ])

# Change column name
colnames(volcano_df)[colnames(volcano_df) == "rel_df[rownames(volcano_df), ]"] <- "rel_abundance"


label_df <- volcano_df %>%
  dplyr::filter(Family %in% unique(plot_abund_df$Family) & Family != "Others" & padj <= 0.05)

# Plot
ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = Significant,
           size = rel_abundance)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(
    data = label_df,
    aes(label = Family),
    size = 3,
    max.overlaps = 50
  ) +
  scale_size_continuous(
    range = c(1, 6),
    name = "relative abundance"
  ) +
  scale_color_manual(values = c(
    "Not Significant" = "grey",
    "Control Enriched" = "#ffac33",
    "Sulfidic Enriched" = "#86b7d3"
  )) +
  theme_bw() +
  labs(
    title = "Volcano plot: Control vs. Sulfidic",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )




################### BRAY CURTIS AND PERMANOVA ON PICRUST2 EC DA RESULTS ################

# bray-curtis on relevant pathways and permanova for significance
relevant_ec_picrust_df_copy <- as.data.frame(relevant_ec_picrust_df_sulf)
# rownames(relevant_ec_picrust_df_copy) <- relevant_ec_picrust_df_copy$description
relevant_ec_picrust_df_copy <- relevant_ec_picrust_df_copy %>% 
  dplyr::select(-c(feature, description))
relevant_ec_picrust_df_copy <- t(relevant_ec_picrust_df_copy)

relevant_rel <- decostand(relevant_ec_picrust_df_copy, method = "total") # tss i.e. rel abundance conversion
bc_dist <- vegdist(relevant_rel, method = "bray")
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  sample_name = rownames(pcoa_res$points),
  PCoA1 = pcoa_res$points[,1],
  PCoA2 = pcoa_res$points[,2]
)
# Calculate variance explained per axis
var_pcoa1 <- round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1)
var_pcoa2 <- round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1)

pcoa_plot_df <- merge(pcoa_df, subset_metadata_sulfidic, by = "sample_name")

# Generate the PCoA Plot
ggplot(pcoa_plot_df, aes(x = PCoA1, y = PCoA2, color = treat)) +
  geom_point(size = 3) +
  labs(x = paste0("PCoA1 (", var_pcoa1, "%)"),
       y = paste0("PCoA2 (", var_pcoa2, "%)"),
       title = "Relevant EC Pathway Abundance Bray-Curtis PCoA") +
  stat_ellipse() +
  theme_minimal()

# Statistical test: PERMANOVA 
# Tests if the centroid/spread of your groups differs significantly
permanova <- adonis2(bc_dist ~ treat, data = subset_metadata_sulfidic, permutations = 999)
print(permanova)

subset_metadata_sulfidic <- as.data.frame(subset_metadata_sulfidic)
rownames(subset_metadata_sulfidic) <- subset_metadata_sulfidic$sample_name
## calculate the betadispersion within each region
betadispers <- betadisper(bc_dist, subset_metadata_sulfidic$treat)

## betadispersion test to see if all regions have the same betadispersion
beta.region = permutest(betadispers)
beta.region

######## SULFIDIC and CONTROL (put into functions for abstraction later) #########
# # bray-curtis on relevant pathways and permanova for significance
# relevant_ec_picrust_df_copy <- as.data.frame(relevant_ec_picrust_df_control)
# # rownames(relevant_ec_picrust_df_copy) <- relevant_ec_picrust_df_copy$description
# relevant_ec_picrust_df_copy <- relevant_ec_picrust_df_copy %>% 
#   dplyr::select(-c(feature, description))
# relevant_ec_picrust_df_copy <- t(relevant_ec_picrust_df_copy)
# 
# 
# relevant_rel <- decostand(relevant_ec_picrust_df_copy, method = "total") # tss i.e. rel abundance conversion
# bc_dist <- vegdist(relevant_rel, method = "bray")
# pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
# pcoa_df <- data.frame(
#   sample_name = rownames(pcoa_res$points),
#   PCoA1 = pcoa_res$points[,1],
#   PCoA2 = pcoa_res$points[,2]
# )
# # Calculate variance explained per axis
# var_pcoa1 <- round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1)
# var_pcoa2 <- round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1)
# 
# pcoa_plot_df <- merge(pcoa_df, subset_metadata_control, by = "sample_name")
# 
# # Generate the PCoA Plot
# ggplot(pcoa_plot_df, aes(x = PCoA1, y = PCoA2, color = treat)) +
#   geom_point(size = 3) +
#   labs(x = paste0("PCoA1 (", var_pcoa1, "%)"),
#        y = paste0("PCoA2 (", var_pcoa2, "%)"),
#        title = "Relevant EC Pathway Abundance Bray-Curtis PCoA") +
#   stat_ellipse() +
#   theme_minimal()
# 
# # Statistical test: PERMANOVA 
# # Tests if the centroid/spread of your groups differs significantly
# permanova <- adonis2(bc_dist ~ treat, data = subset_metadata_control, permutations = 999)
# print(permanova)
# 
# subset_metadata_control <- as.data.frame(subset_metadata_control)
# rownames(subset_metadata_control) <- subset_metadata_control$sample_name
# ## calculate the betadispersion within each region
# betadispers <- betadisper(bc_dist, subset_metadata_control$treat)
# 
# ## betadispersion test to see if all regions have the same betadispersion
# beta.region = permutest(betadispers)
# beta.region


########################## GENE SET ENRICHMENT ANALYSIS ########################

gsea_res <- pathway_gsea(
  abundance = subset_ec,
  metadata = subset_metadata,
  group = "sed_treat",
  pathway_type = "MetaCyc",
  method = "camera",
  seed = 42,
)

## nothing pops up as significant, will try doing with official pipeline in another script
