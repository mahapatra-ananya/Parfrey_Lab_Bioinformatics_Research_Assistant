## load libraries
library(phyloseq)
library(tidyverse)
library(plyr)
library(qualpalr)
library(ggpubr)
library(indicspecies)
library(data.table)
library(DESeq2)
library(ANCOMBC)
library(microbiome)
library(ggrepel)
library(patchwork)

## set wd

# nereocystis_ps_8000 <- readRDS("rarefied_phyloseq_nereocystis_gurkiran_8000reads.RDS")
# nereocystis_ps_20000 <- readRDS("rarefied_phyloseq_nereocystis_gurkiran_20000reads.RDS")
nereocystis_ps <- readRDS("filtered_phyloseq_nereocystis_gurkiran.RDS")
nereocystis_ps_control_WS <- subset_samples(nereocystis_ps, Treatment == "control" & Strain == "WS")
nereocystis_ps_antibiotic_WS <- subset_samples(nereocystis_ps, Treatment == "antibiotic" & Strain == "WS")
nereocystis_ps_control_CR <- subset_samples(nereocystis_ps, Treatment == "control" & Strain == "CR")
nereocystis_ps_antibiotic_CR <- subset_samples(nereocystis_ps, Treatment == "antibiotic" & Strain == "CR")

################################################################
########################### TAXAPLOT ###########################
################################################################

make_combined_df <- function(ps, glom_level, treatment_var) {
  ps_glom = tax_glom(ps, taxrank = glom_level)
  ps_glom@sam_data$read_depth = sample_sums(ps_glom)
  meta = as.data.frame(as.matrix(ps_glom@sam_data)) 
  metacols = ncol(meta)
  otu = as.data.frame(t(as.matrix(ps_glom@otu_table))) |>
    rownames_to_column(var = "label")
  tax = as.data.frame(as.matrix(ps_glom@tax_table)) |>
    rownames_to_column(var = "placeholder_name")
  mo = full_join(meta, otu)
  mo.long = mo |>
    pivot_longer(cols = -c(1:metacols), names_to = "placeholder_name", values_to = "taxa_abundance")
  df = full_join(mo.long, tax, by = "placeholder_name")
  df$relativeabundance = as.numeric(df$taxa_abundance)/as.numeric(df$read_depth)
  if (glom_level == "Family") {
    df$plotnames = df$Family
  } else {
    df$plotnames = paste0(df$Family, ";", df$Genus)
  }
  df.sum = ddply(df, c("plotnames"), summarise, sum = sum(relativeabundance, na.rm = TRUE))
  sorted = df.sum[order(-df.sum$sum), ]
  top = sorted[c(1:15), ]
  top$place = "top"
  print("top:")
  print(top)
  alldata_tops = left_join(top, df)
  print("alldata_tops:")
  print(head(alldata_tops))
  allothers = ddply(alldata_tops, c("label", treatment_var), summarise, top_taxa_sumra = sum(relativeabundance))
  allothers$relativeabundance = 1 - allothers$top_taxa_sumra
  print("allothers:")
  print(head(allothers))
  alldata = full_join(alldata_tops, allothers)
  alldata$place = replace(alldata$place, is.na(alldata$place), "bottom")
  alldata[alldata$place == "bottom", ]$plotnames <- "Others"
  return(alldata)
}

make_taxa_plot <- function(alldata, treatment_var, plotcolors) {
  desired_order <- c("control 10 WS", "control 16 WS", "antibiotic 10 WS", "antibiotic 16 WS",
                     "control 10 CR", "control 16 CR", "antibiotic 10 CR", "antibiotic 16 CR")
  alldata[[treatment_var]] <- factor(alldata[[treatment_var]] , levels = desired_order)
  alldata = alldata[order(-alldata$relativeabundance),]
  natural.order = as.list(c(unique(alldata$plotnames)))
  no.others=natural.order[!natural.order == 'Others']
  plot.order = append(no.others, "Others")
  plot.order = unlist(plot.order)
  alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))
  formula <- paste0(".~", treatment_var)
  p <- ggplot(alldata, aes(x=label, y=as.numeric(relativeabundance),
                      fill=as.factor(plotnames)))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values=plotcolors)+
    guides(fill=guide_legend(ncol=1))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill="white"),
          axis.text.y = element_text(size = 10, colour = "black"),
          axis.title = element_text(size=10, face="bold"),
          strip.text = element_text(color="black", size=10),
          legend.text=element_text(size=6),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 6,
            colour = "black"
          ))+
    labs(y="Relative Abundance", x="Sample", fill="Taxa")+
    facet_grid(formula, scales="free")
  return(p)
}


## Non-rarefied
alldata <- make_combined_df(nereocystis_ps, "Genus", "overall_treat")
all_taxa <- unique(alldata$plotnames)
pal <- qualpal(length(all_taxa), colorspace = "Tableau:20")$hex
plotcolors <- setNames(pal, all_taxa)
plotcolors["Others"] <- "grey90"
taxaplot <- make_taxa_plot(alldata, "overall_treat", plotcolors)
print(taxaplot)
ggsave("taxaplot_nereocystis_gurkiran_nonrare.png", taxaplot)

# ## Rarefied
# ## 20,000 reads
# alldata_20000 <- make_combined_df(nereocystis_ps_20000, "Genus", "overall_treat")
# ## 8,000 reads
# alldata_8000 <- make_combined_df(nereocystis_ps_8000, "Genus", "overall_treat")
# ## Color scheme
# all_taxa <- unique(c(alldata_20000$plotnames, alldata_8000$plotnames))
# set.seed(20)
# pal <- qualpal(length(all_taxa), colorspace = "Tableau:20")$hex
# plotcolors <- setNames(pal, all_taxa)
# plotcolors["Others"] <- "grey90"
# ## Plotting
# taxaplot_20000 <- make_taxa_plot(alldata_20000, "overall_treat", plotcolors)
# print(taxaplot_20000)
# ggsave("taxaplot_nereocystis_gurkiran_20000reads.png", taxaplot_20000)
# taxaplot_8000 <- make_taxa_plot(alldata_8000, "overall_treat", plotcolors)
# print(taxaplot_8000)
# ggsave("taxaplot_nereocystis_gurkiran_8000reads.png", taxaplot_8000)


################################################################
####################### PRESENCE/ABSENCE #######################
################################################################

# Compute indval scores 
compute_indval <- function(ps, comparison_var, comparison_name_for_file) {
  ps@sam_data$read_depth = sample_sums(ps)
  metadata = as.data.frame(ps@sam_data)
  otu = as.data.frame(as.matrix(ps@otu_table))
  otu = as.data.frame(t(as.matrix(ps@otu_table)))
  comparison = metadata[[comparison_var]]
  indval <- multipatt(otu, comparison, duleg = TRUE, control = how(nperm=999))
  indval.str <- as.data.frame(indval$str)
  indval.str$rn <- rownames(indval.str)
  indval.stat <- as.data.frame(indval$sign) 
  indval.stat$rn <- rownames(indval.stat) 
  indval.fid <- as.data.frame(indval$A)
  setDT(indval.fid, keep.rownames = TRUE)[]
  colnames(indval.fid) <- paste0("fid.", colnames(indval.fid))
  names(indval.fid)[names(indval.fid) == 'fid.rn'] <- 'rn'
  indval.prev <- as.data.frame(indval$B)
  setDT(indval.prev, keep.rownames = TRUE)[]
  colnames(indval.prev) <- paste0("prev.", colnames(indval.prev))
  names(indval.prev)[names(indval.prev) == 'prev.rn'] <- 'rn' 
  str.and.stat = full_join(indval.str, indval.stat,
                           by="rn")
  fid.and.spec = full_join(indval.fid, indval.prev,
                           by="rn")
  indval_table = full_join(str.and.stat, fid.and.spec,
                           by="rn")
  tax = as.data.frame(ps@tax_table)
  tax = tax %>% rownames_to_column(var="asv_id")
  names(indval_table)[names(indval_table) == 'rn'] <- 'asv_id'
  indval_table= inner_join(indval_table, tax)
  #write.csv(indval_table, "Results/IndVal/indval_output_", comparison_name_for_file, ".csv")
  metacols= ncol(metadata)+1
  mo = merge(metadata, otu, by=0)
  mo_long = mo |> pivot_longer(cols=-c(1:metacols),
                               names_to = "asv_id",
                               values_to = "asv_abundance")
  sgtaxa = full_join(mo_long, tax)
  sgtaxa$relative_abundance = as.numeric(sgtaxa$asv_abundance)/as.numeric(sgtaxa$read_depth)
  sgtaxa$presabs <- ifelse(sgtaxa$relative_abundance == 0, "asb.", "pres.")
  res <- list(indval_table, sgtaxa)
  return(res)
}

# Make bubble plot
subsetting_for_plot <- function(indval_res, stat_thresh, index_value, taxa_level) {
  indval_table <- indval_res[[1]]
  indval_subset = subset(indval_table,
                         indval_table$stat >= stat_thresh & 
                           indval_table$p.value < 0.05 &
                           indval_table$index == index_value)
  sgtaxa <- indval_res[[2]]
  # Making a variable that indicates if the taxa is enriched in the indval analysis or not
  sgtaxa$enriched = ifelse(sgtaxa$asv_id %in% c(indval_subset$asv_id), "yes", "no")
  
  enriched_plot_df = subset(sgtaxa, sgtaxa$enriched=="yes")
  
  if (taxa_level == "Genus") {
    enriched_plot_df$plotname <- paste0(enriched_plot_df$Family, ";", enriched_plot_df$Genus)
  } else {
    enriched_plot_df$plotname <- enriched_plot_df[[taxa_level]]  
  }
  
  return(enriched_plot_df)
}

make_bubble_plot <- function(enriched_plot_df, stat_thresh, comparison_var_value) {
  # Making bubble plot
  p <- ggplot(enriched_plot_df, aes(x=as.character(Row.names),
                                    y=asv_id,
                                    size=relative_abundance,
                                    alpha=presabs))+
    geom_point()+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 8, colour = "black"),
          axis.text.x = element_blank(),
          axis.title = element_text(size=10),
          strip.text = element_text(color="black", size=8),
          legend.text=element_text(size=8),
          axis.line = element_line(colour = "black"),
          strip.text.y = element_text(angle = 0, size=10))+
    facet_grid(plotname~treat, space="free", scales="free")+
    labs(x="Samples", y = "ASV", size="relative abundance")
  
  #ggsave("Results/IndVal/Indval_EnrichedIn", comparison_var_value, "_stat", stat_thresh, ".png")
}

# not putting below into functions to enable easy inspection of elements 

## Gurkiran Nereocystis:
# ASV level
ind_nereocystis_asv <- compute_indval(nereocystis_ps, "overall_treat", "OverallTreatment_Nereocystis")
head(ind_nereocystis_asv[[1]]) # can also View(ind_rhizome[[1]])
# Genus level
nereocystis_ps_genus <- tax_glom(nereocystis_ps, taxrank = "Genus")
ind_nereocystis_genus <- compute_indval(nereocystis_ps_genus, "overall_treat", "OverallTreatment_Nereocystis")
head(ind_nereocystis_genus[[1]])

# Function to calculate pres abs
make_pres_abs_df <- function(indval_table, pres_col, abs_col, agglom_level) {
  result <- indval_table %>%
    filter({{ pres_col }} > 0 & {{ abs_col }} == 0)
  if (agglom_level == "ASV") {
    result <- result %>%
      mutate(Genus_Species = paste0(Genus, "_", Species)) %>%
      select(asv_id, Genus_Species, {{ pres_col }})
  } else if (agglom_level == "Genus") {
    result <- result %>%
      select(Genus, {{ pres_col }})
  }
  return(result)
}

## ASV level
# pres_10_abs_16_control_WS_df_asv <- make_pres_abs_df(ind_nereocystis_asv[[1]], `prev.control 10 WS`, `prev.control 16 WS`, "ASV")
# pres_10_abs_16_control_CR_df_asv <- make_pres_abs_df(ind_nereocystis_asv[[1]], `prev.control 10 CR`, `prev.control 16 CR`, "ASV")
# pres_16_abs_10_control_WS_df_asv <- make_pres_abs_df(ind_nereocystis_asv[[1]], `prev.control 16 WS`, `prev.control 10 WS`, "ASV")
# pres_16_abs_10_control_CR_df_asv <- make_pres_abs_df(ind_nereocystis_asv[[1]], `prev.control 16 CR`, `prev.control 10 CR`, "ASV")
# pres_10_abs_16_antibiotic_WS_df_asv <- make_pres_abs_df(ind_nereocystis_asv[[1]], `prev.antibiotic 10 WS`, `prev.antibiotic 16 WS`, "ASV")
# pres_10_abs_16_antibiotic_CR_df_asv <- make_pres_abs_df(ind_nereocystis_asv[[1]], `prev.antibiotic 10 CR`, `prev.antibiotic 16 CR`, "ASV")
# pres_16_abs_10_antibiotic_WS_df_asv <- make_pres_abs_df(ind_nereocystis_asv[[1]], `prev.antibiotic 16 WS`, `prev.antibiotic 10 WS`, "ASV")
# pres_16_abs_10_antibiotic_CR_df_asv <- make_pres_abs_df(ind_nereocystis_asv[[1]], `prev.antibiotic 16 CR`, `prev.antibiotic 10 CR`, "ASV")
# 
# combined_prev_df_asv <- ind_nereocystis_asv[[1]] %>%
#   filter(`prev.control 10 WS` == 0 | `prev.control 16 WS` == 0 |
#          `prev.control 10 CR` == 0 | `prev.control 16 CR` == 0|
#          `prev.antibiotic 10 WS` == 0 | `prev.antibiotic 16 WS` == 0|
#          `prev.antibiotic 10 CR` == 0 | `prev.antibiotic 16 CR` == 0) %>%
#   mutate(Genus_Species = paste0(Genus, "_", Species)) %>%
#   select(asv_id, Genus_Species, `prev.control 10 WS`, `prev.control 16 WS`,
#          `prev.control 10 CR`, `prev.control 16 CR`,
#          `prev.antibiotic 10 WS`, `prev.antibiotic 16 WS`,
#          `prev.antibiotic 10 CR`, `prev.antibiotic 16 CR`)
# write.csv(combined_prev_df_asv, "combined_prevalence_df_asv.csv")

# Genus level
pres_10_abs_16_control_WS_df_genus <- make_pres_abs_df(ind_nereocystis_genus[[1]], `prev.control 10 WS`, `prev.control 16 WS`, "Genus")
pres_10_abs_16_control_CR_df_genus <- make_pres_abs_df(ind_nereocystis_genus[[1]], `prev.control 10 CR`, `prev.control 16 CR`, "Genus")
pres_16_abs_10_control_WS_df_genus <- make_pres_abs_df(ind_nereocystis_genus[[1]], `prev.control 16 WS`, `prev.control 10 WS`, "Genus")
pres_16_abs_10_control_CR_df_genus <- make_pres_abs_df(ind_nereocystis_genus[[1]], `prev.control 16 CR`, `prev.control 10 CR`, "Genus")
pres_10_abs_16_antibiotic_WS_df_genus <- make_pres_abs_df(ind_nereocystis_genus[[1]], `prev.antibiotic 10 WS`, `prev.antibiotic 16 WS`, "Genus")
pres_10_abs_16_antibiotic_CR_df_genus <- make_pres_abs_df(ind_nereocystis_genus[[1]], `prev.antibiotic 10 CR`, `prev.antibiotic 16 CR`, "Genus")
pres_16_abs_10_antibiotic_WS_df_genus <- make_pres_abs_df(ind_nereocystis_genus[[1]], `prev.antibiotic 16 WS`, `prev.antibiotic 10 WS`, "Genus")
pres_16_abs_10_antibiotic_CR_df_genus <- make_pres_abs_df(ind_nereocystis_genus[[1]], `prev.antibiotic 16 CR`, `prev.antibiotic 10 CR`, "Genus")

write.csv(pres_10_abs_16_control_WS_df_genus, "pres_10_abs_16_control_WS_genus.csv")
write.csv(pres_10_abs_16_control_CR_df_genus, "pres_10_abs_16_control_CR_genus.csv")
write.csv(pres_16_abs_10_control_WS_df_genus, "pres_16_abs_10_control_WS_genus.csv")
write.csv(pres_16_abs_10_control_CR_df_genus, "pres_16_abs_10_control_CR_genus.csv")
write.csv(pres_10_abs_16_antibiotic_WS_df_genus, "pres_10_abs_16_antibiotic_WS_genus.csv")
write.csv(pres_10_abs_16_antibiotic_CR_df_genus, "pres_10_abs_16_antibiotic_CR_genus.csv")
write.csv(pres_16_abs_10_antibiotic_WS_df_genus, "pres_16_abs_10_antibiotic_WS_genus.csv")
write.csv(pres_16_abs_10_antibiotic_CR_df_genus, "pres_16_abs_10_antibiotic_CR_genus.csv")

combined_prev_df_genus <- ind_nereocystis_genus[[1]] %>%
  filter(`prev.control 10 WS` == 0 | `prev.control 16 WS` == 0 |
         `prev.control 10 CR` == 0 | `prev.control 16 CR` == 0|
         `prev.antibiotic 10 WS` == 0 | `prev.antibiotic 16 WS` == 0|
         `prev.antibiotic 10 CR` == 0 | `prev.antibiotic 16 CR` == 0) %>%
  select(Genus, `prev.control 10 WS`, `prev.control 16 WS`,
         `prev.control 10 CR`, `prev.control 16 CR`,
         `prev.antibiotic 10 WS`, `prev.antibiotic 16 WS`,
         `prev.antibiotic 10 CR`, `prev.antibiotic 16 CR`)
write.csv(combined_prev_df_genus, "combined_prevalence_df_genus.csv")

################################################################
#################### DIFFERENTIAL ABUNDANCE ####################
################################################################

## DESeq2
perform_deseq_analysis <- function(ps, treatment_str) {
  ps_genus <- tax_glom(ps, taxrank = "Genus")
  sample_data(ps_genus)$Temperature <- as.factor(sample_data(ps_genus)$Temperature)
  # Convert ps to DESeq Dataset
  dds = phyloseq_to_deseq2(ps_genus, ~ Temperature)
  # Run the actual differential abundance analysis using Wald test and positive counts (useful for sparse data)
  dds = DESeq(dds, test="Wald", fitType="parametric", sfType="poscounts")
  # Extract results
  # cooksCutoff FALSE so it doesn't filter outliers automatically
  # contrast explicitly defines numerator (R) and denominator (I) for easy interpretation of logfold2 change
  # earlier had it other way around, but switched as intact should be ground truth upon which comparison is based
  res = results(dds, cooksCutoff = FALSE, contrast=c("Temperature", "16", "10"))
  # Table of significant results 
  sigtab = res[which(res$padj < 0.05), ]
  # Merge with taxonomy table (cbind appends columns)
  res = cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  # Save full and significant results
  write.csv(res, paste0("DESeq2_", treatment_str, "_10vs16_FullResults.csv"), row.names = FALSE)
  write.csv(res, paste0("DESeq2_", treatment_str, "_10vs16_SigResults.csv"), row.names = FALSE)
  # Create df for volcano plot
  res <- as.data.frame(res) # convert to df
  volcano_df <- res %>%
    filter(!is.na(padj))
  # Add significance labels to df
  volcano_df$Significant <- "Not Significant"
  volcano_df$Significant[volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 0] <- paste0(treatment_str, " 16C Enriched")
  volcano_df$Significant[volcano_df$padj < 0.05 & volcano_df$log2FoldChange < 0] <- paste0(treatment_str, " 10C Enriched")
  return(volcano_df)
}

# Below is copied from MJ's code for volcano plot from DESeq2 results
# in her 4_microbialanalysis.R file.
# I simply put into a function

create_volcano_plot_deseq2 <- function(ps, volcano_df, treatment_str) {
  ps_genus <- tax_glom(ps, taxrank = "Genus")
  # Calculate relative abundance for each ASV in each sample in OTU table
  ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
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
  #volcano_df$Genus = paste0(volcano_df$Family, ";", volcano_df$Genus)
  label_df <- volcano_df %>%
    filter(padj < 0.05)
  colors <- c(
    "Not Significant" = "grey",
    "Temp1" = "red",
    "Temp2" = "blue"
  )
  # Dynamically rename the keys using your 'treatment' variable
  names(colors)[2] <- paste0(treatment_str, " 16C Enriched")
  names(colors)[3] <- paste0(treatment_str, " 10C Enriched")
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
      size = 2.5,
      max.overlaps = Inf
    ) +
    scale_size_continuous(
      range = c(1, 6),
      name = "relative abundance"
    ) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      title = paste0("Volcano plot: ", treatment_str, " 10C vs 16C"),
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    )
  ggsave(paste0("DESeq2_Volcano_", treatment_str, "_10vs16.png"), plot = p)
  return(p)
}

# Call functions and look at plots
## Agglomerated at Genus Level 
#Control WS 
controlWS_deseq2_results_genus <- perform_deseq_analysis(nereocystis_ps_control_WS, "Control_WS")
controlWS_deseq2_volcanoplot_genus <- create_volcano_plot_deseq2(nereocystis_ps_control_WS, controlWS_deseq2_results_genus, "Control_WS")
controlWS_deseq2_volcanoplot_genus 
# Control CR 
controlCR_deseq2_results_genus <- perform_deseq_analysis(nereocystis_ps_control_CR, "Control_CR")
controlCR_deseq2_volcanoplot_genus <- create_volcano_plot_deseq2(nereocystis_ps_control_CR, controlCR_deseq2_results_genus, "Control_CR")
controlCR_deseq2_volcanoplot_genus 
# Antibiotic WS
antibioticWS_deseq2_results_genus <- perform_deseq_analysis(nereocystis_ps_antibiotic_WS, "Antibiotic_WS")
antibioticWS_deseq2_volcanoplot_genus <- create_volcano_plot_deseq2(nereocystis_ps_antibiotic_WS, antibioticWS_deseq2_results_genus, "Antibiotic_WS")
antibioticWS_deseq2_volcanoplot_genus 
# Antibiotic CR
antibioticCR_deseq2_results_genus <- perform_deseq_analysis(nereocystis_ps_antibiotic_CR, "Antibiotic_CR")
antibioticCR_deseq2_volcanoplot_genus <- create_volcano_plot_deseq2(nereocystis_ps_antibiotic_CR, antibioticCR_deseq2_results_genus, "Antibiotic_CR")
antibioticCR_deseq2_volcanoplot_genus

combined_deseq2_plot <- ((controlWS_deseq2_volcanoplot_genus + antibioticWS_deseq2_volcanoplot_genus) /
                            (controlCR_deseq2_volcanoplot_genus + antibioticCR_deseq2_volcanoplot)) + plot_layout(axes = "collect_x")
combined_deseq2_plot
ggsave("deseq2_combined_volcanoplot_genus.png", combined_deseq2_plot)

## ASV level
# # Control WS
# controlWS_deseq2_results <- perform_deseq_analysis(nereocystis_ps_control_WS, "Control_WS")
# controlWS_deseq2_volcanoplot <- create_volcano_plot_deseq2(nereocystis_ps_control_WS, controlWS_deseq2_results, "Control_WS")
# controlWS_deseq2_volcanoplot
# # Control CR
# controlCR_deseq2_results <- perform_deseq_analysis(nereocystis_ps_control_CR, "Control_CR")
# controlCR_deseq2_volcanoplot <- create_volcano_plot_deseq2(nereocystis_ps_control_CR, controlCR_deseq2_results, "Control_CR")
# controlCR_deseq2_volcanoplot
# # Antibiotic WS
# antibioticWS_deseq2_results <- perform_deseq_analysis(nereocystis_ps_antibiotic_WS, "Antibiotic_WS")
# antibioticWS_deseq2_volcanoplot <- create_volcano_plot_deseq2(nereocystis_ps_antibiotic_WS, antibioticWS_deseq2_results, "Antibiotic_WS")
# antibioticWS_deseq2_volcanoplot
# # Antibiotic CR
# antibioticCR_deseq2_results <- perform_deseq_analysis(nereocystis_ps_antibiotic_CR, "Antibiotic_CR")
# antibioticCR_deseq2_volcanoplot <- create_volcano_plot_deseq2(nereocystis_ps_antibiotic_CR, antibioticCR_deseq2_results, "Antibiotic_CR")
# antibioticCR_deseq2_volcanoplot

## ANCOM-BC
get_ancombc_results <- function(ps) {
  sample_data(ps)$Temperature <- as.factor(sample_data(ps)$Temperature)
  out = ancombc(data = ps, tax_level = "Genus", p_adj_method = "BH",
                    formula = "Temperature", group = "Temperature")
  res = out$res
  return(res)
}

create_ancombc_volcano_df <- function(ps, res) {
  ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
  ps_rel   <- transform_sample_counts(ps_genus, function(x) x / sum(x))
  rel_mat <- as.data.frame(as(otu_table(ps_rel), "matrix"))
  if (!taxa_are_rows(ps_rel)) {
    rel_mat <- t(rel_mat)
  }
  tax_mat <- as.data.frame(as(tax_table(ps_rel), "matrix"))
  rel_df <- data.frame(
    Genus = tax_mat$Genus,
    rel_abundance = rowMeans(rel_mat),
    stringsAsFactors = FALSE
  )
  rel_df <- aggregate(rel_abundance ~ Genus, data = rel_df, FUN = mean)
  volcano_df <- data.frame(
    Genus = res$lfc$taxon,
    LFC = res$lfc$Temperature16,
    p_val = res$p_val$Temperature16,
    q_val = res$q_val$Temperature16,
    stringsAsFactors = FALSE
  )
  volcano_df <- volcano_df[!is.na(volcano_df$q_val), ]
  volcano_df$Significant <- "Not Significant"
  volcano_df$Significant[volcano_df$q_val < 0.05 & volcano_df$LFC > 0] <- "16C Enriched"
  volcano_df$Significant[volcano_df$q_val < 0.05 & volcano_df$LFC < 0] <- "10C Enriched"
  volcano_df$Significant <- factor(
    volcano_df$Significant,
    levels = c("Not Significant", "16C Enriched", "10C Enriched")
  )
  volcano_df <- merge(volcano_df, rel_df, by = "Genus", all.x = TRUE)
  return(volcano_df)
}

create_volcanoplot_ancombc <- function(volcano_df, treatment_str) {
  label_df <- volcano_df[volcano_df$q_val < 0.05, ]
  label_df <- label_df[!duplicated(label_df$Genus), ]
  color_map <- c(
    "Not Significant" = "grey",
    "16C Enriched" = "red",
    "10C Enriched" = "blue"
  )
  max_neg_log10p <- max(-log10(volcano_df$q_val), na.rm = TRUE)
  max_abs_lfc <- max(abs(volcano_df$LFC), na.rm = TRUE)
  p <- ggplot(volcano_df, aes(
    x = LFC,
    y = -log10(q_val),
    color = Significant,
    size = rel_abundance
  )) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(range = c(1, 6), name = "Relative Abundance") +
    geom_text_repel(
      data = label_df,
      aes(label = Genus),
      size = 2.5,
      max.overlaps = Inf,
      show.legend = FALSE,
      box.padding = 0.4,
      point.padding = 0.2,
      segment.color = "grey50",
      segment.size = 0.2,
      force = 10
    ) +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_x_continuous(
      limits = c(-max_abs_lfc, max_abs_lfc),
      expand = expansion(mult = 0.15)
    ) +
    scale_y_continuous(
      limits = c(0, max_neg_log10p),
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_bw() +
    labs(
      title = paste0(treatment_str, ": ANCOMBC 10C vs. 16C"),
      x = "Log-Fold Change (Natural Log scale)",
      y = "-log10 p-adjusted value"
    )
  ggsave(paste0("ANCOMBC_VolcanoPlot_", treatment_str, "_10Cvs16C_Genus.png"), p)
  return(p)
}

# Control WS
ancombc_res_control_WS <- get_ancombc_results(nereocystis_ps_control_WS)
ancombc_volcanodf_control_WS <- create_ancombc_volcano_df(nereocystis_ps_control_WS, ancombc_res_control_WS)
ancombc_volcanoplot_control_WS <- create_volcanoplot_ancombc(ancombc_volcanodf_control_WS, "Control_WS")
ancombc_volcanoplot_control_WS

# Control CR
ancombc_res_control_CR <- get_ancombc_results(nereocystis_ps_control_CR)
ancombc_volcanodf_control_CR <- create_ancombc_volcano_df(nereocystis_ps_control_CR, ancombc_res_control_CR)
ancombc_volcanoplot_control_CR <- create_volcanoplot_ancombc(ancombc_volcanodf_control_CR, "Control_CR")
ancombc_volcanoplot_control_CR

# Antibiotic WS
ancombc_res_antibiotic_WS <- get_ancombc_results(nereocystis_ps_antibiotic_WS)
ancombc_volcanodf_antibiotic_WS <- create_ancombc_volcano_df(nereocystis_ps_antibiotic_WS, ancombc_res_antibiotic_WS)
ancombc_volcanoplot_antibiotic_WS <- create_volcanoplot_ancombc(ancombc_volcanodf_antibiotic_WS, "Antibiotic_WS")
ancombc_volcanoplot_antibiotic_WS

# Antibiotic CR
ancombc_res_antibiotic_CR <- get_ancombc_results(nereocystis_ps_antibiotic_CR)
ancombc_volcanodf_antibiotic_CR <- create_ancombc_volcano_df(nereocystis_ps_antibiotic_CR, ancombc_res_antibiotic_CR)
ancombc_volcanoplot_antibiotic_CR <- create_volcanoplot_ancombc(ancombc_volcanodf_antibiotic_CR, "Antibiotic_CR")
ancombc_volcanoplot_antibiotic_CR

combined_ancombc_plot <- ((ancombc_volcanoplot_control_WS + ancombc_volcanoplot_antibiotic_WS) /
  (ancombc_volcanoplot_control_CR + ancombc_volcanoplot_antibiotic_CR)) + plot_layout(axes = "collect_x")
combined_ancombc_plot
ggsave("ancombc_combined_volcanoplot_genus.png", combined_ancombc_plot)
