# Load libraries
library(phyloseq)
library(tidyverse)
library(indicspecies)
library(data.table)

## set wd

# Read in data and subset
seagrass <- readRDS("../../Data/filtered_seagrass_mesocosm.RDS")
seagrass.rhizome.day20 <- subset_samples(seagrass, sample_type == "rhizome" & experiment_day == 20)
seagrass.control <- subset_samples(seagrass.rhizome.day20, sed_treat == "control")
seagrass.sulfidic <- subset_samples(seagrass.rhizome.day20, sed_treat == "sulfidic")

# Compute indval scores and make bubble plot
compute_indval <- function(ps, comparison_var, comparison_name_for_file) {
  # Calculate sample sums (read depth) of dataset
  ps@sam_data$read_depth = sample_sums(ps)
  
  # Extracting data tables
  metadata = as.data.frame(ps@sam_data)
  otu = as.data.frame(as.matrix(ps@otu_table))
  # Have to transpose otu table to get samples as rows and ASVs as columns
  otu = as.data.frame(t(as.matrix(ps@otu_table)))
  
  # Creating an object with the comparison stored
  comparison = metadata[[comparison_var]]
  
  # Running indval 
  indval <- multipatt(otu, comparison, duleg = TRUE, control = how(nperm=999))
  
  # Extracting indval statistic from indval object
  indval.str <- as.data.frame(indval$str)
  indval.str$rn <- rownames(indval.str)
  
  # Extracting p-value and index value
  indval.stat <- as.data.frame(indval$sign) #get dataframe of indval statistic
  indval.stat$rn <- rownames(indval.stat) # make column of ASVs
  
  # Extracting Fidelity as dataframe
  indval.fid <- as.data.frame(indval$A)
  # Extracting rownames into column
  setDT(indval.fid, keep.rownames = TRUE)[]
  
  # Renaming columns
  colnames(indval.fid) <- paste0("fid.", colnames(indval.fid))
  names(indval.fid)[names(indval.fid) == 'fid.rn'] <- 'rn'
  
  # Extracting Prevalence as dataframe
  indval.prev <- as.data.frame(indval$B)
  # Extracting rownames into column
  setDT(indval.prev, keep.rownames = TRUE)[]
  
  # Renaming columns
  colnames(indval.prev) <- paste0("prev.", colnames(indval.prev))
  names(indval.prev)[names(indval.prev) == 'prev.rn'] <- 'rn'
  
  # Joining statistics together using the rn column found in each dataframe. 
  # We use the rn column for merging the stats, fidelity, and prevalence output 
  # into one table because this column is identical between every dataframe. 
  str.and.stat = full_join(indval.str, indval.stat,
                           by="rn")
  fid.and.spec = full_join(indval.fid, indval.prev,
                           by="rn")
  indval_table = full_join(str.and.stat, fid.and.spec,
                           by="rn")
  
  
  # Extracting the taxonomy table from the phyloseq object
  tax = as.data.frame(ps@tax_table)
  # Creating a new column containing the ASV IDs, which I will use to merge the dataframes.
  tax = tax %>% rownames_to_column(var="asv_id")
  
  # Renaming columns from the IndVal output to join with taxonomy
  names(indval_table)[names(indval_table) == 'rn'] <- 'asv_id'
  
  # Merging with taxonomy
  indval_table= inner_join(indval_table, tax)
  
  # Save table
  #write.csv(indval_table, "Results/IndVal/indval_output_", comparison_name_for_file, ".csv")
  
  metacols= ncol(metadata)+1
  # Combining otu and metadata
  mo = merge(metadata, otu, by=0)
  # Pivoting table to add taxonomy
  mo_long = mo |> pivot_longer(cols=-c(1:metacols),
                               names_to = "asv_id",
                               values_to = "asv_abundance")
  # Joining with taxonomy
  sgtaxa = full_join(mo_long, tax)
  
  # Calculating the relative abundance of each ASV in each sample
  sgtaxa$relative_abundance = as.numeric(sgtaxa$asv_abundance)/as.numeric(sgtaxa$read_depth)
  # Making a presence/absence variable for each asv in each sample
  sgtaxa$presabs <- ifelse(sgtaxa$relative_abundance == 0, "asb.", "pres.")
  
  res <- list(indval_table, sgtaxa)
  return(res)
}

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

# Whole rhizome:
# compute indval table
ind_rhizome <- compute_indval(seagrass, "sample_type", "RhizomevsBulk")
# view full indval table
head(ind_rhizome[[1]]) # can also View(ind_rhizome[[1]])
# pass full res into subset for plot function, look for enriched in rhizome
enriched_df_rhizome <- subsetting_for_plot(ind_rhizome, 0.8, 2, "Family")
# subset table to only rhizome samples
enriched_df_rhizome <- subset(enriched_df_rhizome, enriched_df_rhizome$sample_type=="rhizome" & experiment_day==20)
# pass enriched df into plotting function
plot_rhizome <- make_bubble_plot(enriched_df_rhizome, 0.8, "RhizomeSamples")
# View plot
plot_rhizome

# Sediment treatment: 
# compute indval table
ind_sediment <- compute_indval(seagrass.rhizome.day20, "sed_treat", "SvsC")
# view full indval table
head(ind_sediment[[1]]) # can also View(ind_sediment[[1]])
# pass full res into subset for plot function, look for enriched in control
enriched_df_control <- subsetting_for_plot(ind_sediment, 0.9, 1, "Family")
# subset table to only control sediment
enriched_df_control <- subset(enriched_df_control, enriched_df_control$sed_treat=="control")
# pass enriched df into plotting function
plot_control <- make_bubble_plot(enriched_df_control, 0.9, "ControlSediment")
# View plot
plot_control
# pass full res into subset for plot function, look for enriched in sulfidic
enriched_df_sulfidic <- subsetting_for_plot(ind_sediment, 0.95, 2, "Family")
# subset to only sulfidic sediment
enriched_df_sulfidic <- subset(enriched_df_sulfidic, enriched_df_sulfidic$sed_treat=="sulfidic")
# pass enriched df into plotting function
plot_sulfidic <- make_bubble_plot(enriched_df_sulfidic, 0.95, "SulfidicSediment")
# View plot
plot_sulfidic

# Sulfidic microbe treatment:
ind_sulfidic_microbe <- compute_indval(seagrass.sulfidic, "microbe_treat", "SIvsSR")
# view full indval table
head(ind_sulfidic_microbe[[1]]) # can also View(ind_sulfidic_microbe[[1]])
# pass full res into subset for plot function, look for enriched in intact
enriched_df_SI <- subsetting_for_plot(ind_sulfidic_microbe, 0.7, 1, "Family")
# subset table to only intact
enriched_df_SI <- subset(enriched_df_SI, enriched_df_SI$microbe_treat=="intact")
# pass enriched df into plotting function
plot_SI <- make_bubble_plot(enriched_df_SI, 0.7, "SulfidicIntact")
# View plot
plot_SI
# pass full res into subset for plot function, look for enriched in removed
enriched_df_SR <- subsetting_for_plot(ind_sulfidic_microbe, 0.7, 2, "Family")
# subset table to only removed
enriched_df_SR <- subset(enriched_df_SR, enriched_df_SR$microbe_treat=="removed")
# pass enriched df into plotting function
plot_SR <- make_bubble_plot(enriched_df_SR, 0.7, "SulfidicRemoved")
# View plot
plot_SR

# Contro microbe treatment
ind_control_microbe <- compute_indval(seagrass.control, "microbe_treat", "CIvsCR")
# view full indval table
head(ind_control_microbe[[1]]) # can also View(ind_control_microbe[[1]])
# pass full res into subset for plot function, look for enriched in intact
enriched_df_CI <- subsetting_for_plot(ind_control_microbe, 0.7, 1, "Family")
# subset table to only intact
enriched_df_CI <- subset(enriched_df_CI, enriched_df_CI$microbe_treat=="intact")
# pass enriched df into plotting function
plot_CI <- make_bubble_plot(enriched_df_CI, 0.7, "ControlIntact")
# View plot
plot_CI
# pass full res into subset for plot function, look for enriched in removed
enriched_df_CR <- subsetting_for_plot(ind_control_microbe, 0.7, 2, "Family")
# subset table to only removed
enriched_df_CR <- subset(enriched_df_CR, enriched_df_CR$microbe_treat=="removed")
# pass enriched df into plotting function
plot_CR <- make_bubble_plot(enriched_df_CR, 0.7, "ControlRemoved")
# View plot
plot_CR
