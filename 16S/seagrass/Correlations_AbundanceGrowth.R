## load libraries
library(phyloseq)
library(tidyverse)
library(vegan)
library(car)
library(ggthemes)
library(plyr)
library(ggeffects)
library(ggpubr)

## set working directory 

# GLOBAL VARIABLES
# Read in data
shootgrowth <- read.csv("../../Data/seagrass_L3growth.csv")
biomass <- read.csv("../../Data/sg_biomass.csv")
# seagrass = readRDS("../../Data/rarefied_seagrass_mesocosm.RDS")
seagrass = readRDS("../../Data/filtered_seagrass_mesocosm.RDS")

# HELPER FUNCTIONS

# Function to agglomerate dataset by level specified
agglomerate <- function(level) {
  copy = tax_glom(seagrass, taxrank = level)
  return(copy)
}

#filter_prevalence <- function(ps, threshold = 0.5) {
 # otu <- as(otu_table(ps), "matrix")
  ## ensure taxa are rows
  #if (!taxa_are_rows(ps)) {
   # otu <- t(otu)
  #}
  #prevalence <- rowSums(otu > 0) / ncol(otu)
  #keep_taxa <- names(prevalence[prevalence >= threshold])
  #prune_taxa(keep_taxa, ps)
#}

# Function to extract phyloseq tables and merge all together
merge_data <- function (agglomerated_phyloseq) {
  # calculate number of reads per sample and put in read depth column
  agglomerated_phyloseq@sam_data$read_depth = sample_sums(agglomerated_phyloseq)
  # extract metadata
  meta = as.data.frame(as.matrix(agglomerated_phyloseq@sam_data)) |>
    rownames_to_column(var = "sample_id")
  metacols = ncol(meta)
  # extract otu
  otu = as.data.frame(t(as.matrix(agglomerated_phyloseq@otu_table))) |>
    rownames_to_column(var = "sample_name")
  # extract taxonomy
  tax = as.data.frame(as.matrix(agglomerated_phyloseq@tax_table)) |>
    rownames_to_column(var = "placeholder_name")
  ## combine metadata and otu data
  mo = full_join(meta, otu, by = "sample_name")
  ## pivot mo longer to be able to join with taxonomy
  mo.long = mo |>
    pivot_longer(cols = -c(1:metacols), names_to = "placeholder_name", values_to = "taxa_abundance")
  ## combine with taxonomy
  mesocosmdf = full_join(mo.long, tax)
  ## calculate the relative abundance of each taxa at agglomerated level within each sample
  mesocosmdf$relativeabundance = as.numeric(mesocosmdf$taxa_abundance)/as.numeric(mesocosmdf$read_depth)
  # renaming column for join
  colnames(mesocosmdf)[colnames(mesocosmdf) == "core_id"] <- "shoot_id"
  # merging full merged metadata, otu and taxa table with shootgrowth and biomass data
  merged_df <- merge(mesocosmdf, shootgrowth, by = "shoot_id", all = FALSE)
  merged_df <- merge(merged_df, biomass, by = "shoot_id", all = FALSE)
  return(merged_df)
}

# Function to iterate through all taxa at specified level and compute correlations
make_final_corr_df <- function(df, level, var) {
  # initialize result, iterate through all taxa and compute correlations by subsetting
  result <- c()
  for (taxa in unique(df[[level]])) {
    corr_list <- all_corr(df, level, taxa, var)
    result <- rbind(result, corr_list)
  }
  # create interpretable column names
  colnames(result) <- c(level, "CI r", "CI p", "CR r", "CR p",
                                   "SI r", "SI p", "SR r", "SR p",
                                   "Overall r", "Overall p")
  return(as.data.frame(result, stringsAsFactors = FALSE))
}

# Function to compute prevalence within the treatment
compute_prevalence <- function(df) {
  # sums the TRUE values (rows where taxa_abundance>0) and divides by total samples
  sum(df$taxa_abundance > 0, na.rm = TRUE) / length(df$taxa_abundance)
  
}

all_corr <- function(input_df, level, taxa, var) {
  # initialize result and obtain subset dfs for each treatment and for the taxa
  result <- c(taxa)
  subsets <- perform_subset(input_df, level, taxa)
  # compute prevalence then compute correlations for all 5 dfs and append to result
  for (df in subsets) {
    prevalence <- compute_prevalence(df)
    if (prevalence < 0.5) {
      result <- c(result, NA, NA) # don't compute corr if below 0.5 prevalence in the treatment
    } else {
      corr <- compute_corr(df$relativeabundance, df[[var]])
      result <- c(result, corr[1], corr[2])
    }
  }
  return(result)
}

# Function to create list of subset dfs by treatments, for a specified taxa
perform_subset <- function(input_df, level, taxa) {
  # subset by taxa
  subset_df = input_df[input_df[[level]] == taxa, ]
  # subset further, once per treatment
  si = subset(subset_df, treat.x == "sulfidic_intact")
  sr = subset(subset_df, treat.x == "sulfidic_removed")
  ci = subset(subset_df, treat.x == "control_intact")
  cr = subset(subset_df, treat.x == "control_removed")
  # combine into list
  result <- list(ci, cr, si, sr, subset_df)
  return(result)
}

# Function to compute correlation between two variables and output r and p-value
compute_corr <- function(v1, v2) {
  v1 <- as.numeric(v1)
  v2 <- as.numeric(v2)
  corr <- cor.test(v1, v2)
  r <- corr$estimate
  p <- corr$p.value
  result <- c(r, p)
  return(result)
}

# MAIN

main <- function(level, var) {
  # agglomerate and merge
  agglom <- agglomerate(level)
  merged_df <- merge_data(agglom)
  # compute correlations
  result <- make_final_corr_df(merged_df, level, var)
  return(result)
}

# Run functions (could reduce computing time by returning all results in a list
# and indexing, to avoid calling functions multiple times - later TODO)
RGR_corr_df_family <- main("Family", "L3_RGR")
biomass_corr_df_family <- main("Family", "rhizome_biomass")
RGR_corr_df_genus <- main("Genus", "L3_RGR")
biomass_corr_df_genus <- main("Genus", "rhizome_biomass")

# Save outputs
# write.csv(RGR_corr_df_genus, "Results/Correlations/Rarefied/ShootGrowth_Abundance_Correlations_ByGenus.csv", row.names = FALSE)
# write.csv(biomass_corr_df_genus, "Results/Correlations/Rarefied/RhizomeBiomass_Abundance_Correlations_ByGenus.csv", row.names = FALSE)
# write.csv(RGR_corr_df_family, "Results/Correlations/Rarefied/ShootGrowth_Abundance_Correlations_ByFamily.csv", row.names = FALSE)
# write.csv(biomass_corr_df_family, "Results/Correlations/Rarefied/RhizomeBiomass_Abundance_Correlations_ByFamily.csv", row.names = FALSE)

write.csv(RGR_corr_df_genus, "Results/Correlations/NonRarefied/PrevFiltered/ShootGrowth_Abundance_Correlations_ByGenus.csv", row.names = FALSE)
write.csv(biomass_corr_df_genus, "Results/Correlations/NonRarefied/PrevFiltered/RhizomeBiomass_Abundance_Correlations_ByGenus.csv", row.names = FALSE)
write.csv(RGR_corr_df_family, "Results/Correlations/NonRarefied/PrevFiltered/ShootGrowth_Abundance_Correlations_ByFamily.csv", row.names = FALSE)
write.csv(biomass_corr_df_family, "Results/Correlations/NonRarefied/PrevFiltered/RhizomeBiomass_Abundance_Correlations_ByFamily.csv", row.names = FALSE)



############## CORRELATION PLOTS ############

make_plots <- function(input_df, taxa, level, feature, sediment) {
  
  subset_df <- input_df[input_df[[level]] == taxa, ]
  subset_df <- subset_df[subset_df$sed_treat.x == sediment, ]
  subset_df$presence <- subset_df$relativeabundance>0 # outputs 1 or 0
  
  results_list <- list()
  
  corr_plot <- ggplot(subset_df,
         aes(x = as.numeric(relativeabundance),
             y = as.numeric(.data[[feature]]),
             color = treat.x)) +
    geom_point(size = 5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 3) +
    theme_classic() +
    labs(
      x = paste("Relative Abundance of Rhizome", taxa),
      y = paste(feature),
      color = "Treatment"
    ) +
    theme(
      axis.title = element_text(size = 27),
      legend.position = "right",
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 30)
    )
  
  box_plot <- ggplot(subset_df,
                     aes(x = treat.x,
                         y = as.numeric(.data$relativeabundance),
                         fill = treat.x)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width = 0.2, height = 0) +
    #geom_point() +
    stat_compare_means(method = "wilcox.test") +
    theme_classic() +
    labs(x = "Treatment",
         y = paste("Relative Abundance of Rhizome", taxa))
  
  test <- fisher.test(table(subset_df$treat.x, subset_df$presence))
  prev_SI <- mean(subset_df$presence[subset_df$treat.x == "sulfidic_intact"])
  prev_SR <- mean(subset_df$presence[subset_df$treat.x == "sulfidic_removed"])
  
  results_list[["correlation"]] = corr_plot
  results_list[["abundance_box"]] = box_plot
  results_list[["presence_fisher"]] = test
  results_list[["prevalence"]] = list("SI" = prev_SI, "SR" = prev_SR)
  
  return(results_list)
}

ahrensia_plots <- make_plots(biomass_corr_df_genus[[1]], "Ahrensia", "Genus", "rhizome_biomass", "sulfidic")
syntrophomonas_plots <- make_plots(biomass_corr_df_genus[[1]], "Syntrophomonas", "Genus", "rhizome_biomass", "sulfidic")
ferrimonas_plots <- make_plots(biomass_corr_df_genus[[1]], "Ferrimonas", "Genus", "rhizome_biomass", "sulfidic")
ahrensia_plots[["abundance_box"]]
syntrophomonas_plots[["abundance_box"]]
ferrimonas_plots[["abundance_box"]]

ahrensia_plots[["presence_fisher"]] # p-value 1; 95% CI 0.0067986 10.5137419
syntrophomonas_plots[["presence_fisher"]] # p-value 0.3498; 95% CI 0.4048219 53.8102540
ferrimonas_plots[["presence_fisher"]] # p-value 0.6285; 95% CI 0.02689425 3.88582019

ahrensia_plots[["prevalence"]]
syntrophomonas_plots[["prevalence"]]
ferrimonas_plots[["prevalence"]]

prevalence_df <- data.frame(
  SI_Prevalence = c(ahrensia_plots[["prevalence"]]$SI, syntrophomonas_plots[["prevalence"]]$SI, ferrimonas_plots[["prevalence"]]$SI),
  SR_Prevalence = c(ahrensia_plots[["prevalence"]]$SR, syntrophomonas_plots[["prevalence"]]$SI, ferrimonas_plots[["prevalence"]]$SR),
  row.names = c("Ahrensia", "Syntrophomonas", "Ferrimonas")
)

prevalence_df
