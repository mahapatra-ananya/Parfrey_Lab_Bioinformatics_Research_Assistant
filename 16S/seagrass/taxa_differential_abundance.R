library(phyloseq)
library("DESeq2")
library(ggrepel)
library(ggplot2)
library(dplyr)

## set wd

# Read in data and subset
seagrass <- readRDS("../../Data/filtered_seagrass_mesocosm.RDS")
seagrass.subset.sulfidic <- subset_samples(seagrass, sample_type == "rhizome" & experiment_day == 20 & sed_treat == "sulfidic")
seagrass.subset.control <- subset_samples(seagrass, sample_type == "rhizome" & experiment_day == 20 & sed_treat == "control")

# Implement functions to perform deseq2 analysis and create plots (saving tables and plots into results folder)
perform_deseq_analysis <- function(ps, sediment) {
  # Convert ps to DESeq Dataset
  dds = phyloseq_to_deseq2(ps, ~ microbe_treat)
  # Run the actual differential abundance analysis using Wald test and positive counts (useful for sparse data)
  dds = DESeq(dds, test="Wald", fitType="parametric", sfType="poscounts")
  # Extract results
  # cooksCutoff FALSE so it doesn't filter outliers automatically
  # contrast explicitly defines numerator (R) and denominator (I) for easy interpretation of logfold2 change
  # earlier had it other way around, but switched as intact should be ground truth upon which comparison is based
  res = results(dds, cooksCutoff = FALSE, contrast=c("microbe_treat", "removed", "intact"))
  # Table of significant results 
  sigtab = res[which(res$padj < 0.05), ]
  # Merge with taxonomy table (cbind appends columns)
  res = cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  # Save full and significant results
  write.csv(res, paste0("Results/DESeq2_", sediment, "_IvsR_FullResults.csv"), row.names = FALSE)
  write.csv(res, paste0("Results/DESeq2_", sediment, "_IvsR_SigResults.csv"), row.names = FALSE)
  # Create df for volcano plot
  res <- as.data.frame(res) # convert to df
  volcano_df <- res %>%
    filter(!is.na(padj))
  # Add significance labels to df
  volcano_df$Significant <- "Not Significant"
  volcano_df$Significant[volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 0] <- paste0(sediment, " Removed Enriched")
  volcano_df$Significant[volcano_df$padj < 0.05 & volcano_df$log2FoldChange < 0] <- paste0(sediment, " Intact Enriched")
  return(volcano_df)
}


# Below is copied from MJ's code for volcano plot from DESeq2 results
# I simply put into a function
# In her 4_microbialanalysis.R file

create_volcano_plot <- function(ps, volcano_df, sediment) {
  # Calculate relative abundance for each ASV in each sample in OTU table
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
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
  volcano_df$Genus = paste0(volcano_df$Family, ";", volcano_df$Genus)
  label_df <- volcano_df %>%
    filter(padj < 0.05)
  colors <- c(
    "Not Significant" = "grey",
    "Temp1" = "red",
    "Temp2" = "blue"
  )
  # Dynamically rename the keys using your 'sediment' variable
  names(colors)[2] <- paste0(sediment, " Removed Enriched")
  names(colors)[3] <- paste0(sediment, " Intact Enriched")
  print(colors)
  # Plot
  p <- ggplot(volcano_df,
         aes(x = log2FoldChange,
             y = -log10(padj),
             color = Significant,
             size = rel_abundance)) +
    geom_point(alpha = 0.8) +
    geom_text_repel(
      data = label_df,
      aes(label = Genus),
      size = 3,
      max.overlaps = 50
    ) +
    scale_size_continuous(
      range = c(1, 6),
      name = "relative abundance"
    ) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      title = paste0("Volcano plot: ", sediment, " Intact vs ", sediment, " Removed"),
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    )
  ggsave(paste0("Results/DESeq2_Volcano_", sediment, "_IvsR.png"), plot = p)
  return(p)
}

# Call functions and look at plots
sulfidic_deseq2_results <- perform_deseq_analysis(seagrass.subset.sulfidic, "Sulfidic")
control_deseq2_results <- perform_deseq_analysis(seagrass.subset.control, "Control")
sulfidic_volcano_plot <- create_volcano_plot(seagrass.subset.sulfidic, sulfidic_deseq2_results, "Sulfidic")
control_volcano_plot <- create_volcano_plot(seagrass.subset.control, control_deseq2_results, "Control")
sulfidic_volcano_plot
control_volcano_plot
