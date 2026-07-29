library(phyloseq)
library(dplyr)
library(stringr)

## set wd

ps <- readRDS("raw_phyloseq_gametophyte.RDS")

## Remove degree symbol in experiment as R does not read it
# meta <- sample_data(ps)
# # meta$experiment <- iconv(meta$experiment, from = "latin1", to = "UTF-8", sub = "byte")
# # #replacements <- c("")
# # meta$experiment <- sub("(16|10).C", "\\1C", meta$experiment)
# #sub("10.C", "10C", meta$experiment)
# sample_data(ps) <- meta
# View(sample_data(ps))
# saveRDS(ps, "raw_phyloseq_character_removed.RDS")

################################################################## 
############################ FILTERING ########################### 
################################################################## 

## 1. Remove off target taxa and contaminant
ps = subset_taxa(ps, Kingdom != "Unassigned" & Order != "Chloroplast" &
                        Family != "Mitochondria" & Kingdom != "Eukaryota" & Genus!= 	
                        "Escherichia-Shigella")
View(ps@tax_table)

## 2. Remove samples with low total reads
sample_sums(ps)
plot(sort(sample_sums(ps)))
ps@sam_data$sample_sums_unfiltered = as.numeric(sample_sums(ps))
View(sample_data(ps))
## All reads >1000, so keeping all

## 3. Remove ASVs found less than 100 times and ASV contaminants
## Less than 100 times in whole dataset
otutab <- as.data.frame(t(as.matrix(otu_table(ps@otu_table))))
## ASVs should be ROWS an sample names COLUMNS. 
View(otutab) # Yes
otutab$asv_abundance = rowSums(otutab)
min(otutab$asv_abundance) # 0
otu.pruned = subset(otutab, otutab$asv_abundance >= 100)
min(otu.pruned$asv_abundance) # Now at 100
widthotu = ncol(otu.pruned)
otu.pruned = otu.pruned[, -c(widthotu)]
## Contaminants (less than 3 times in a sample)
ASVoccur = function(x) {
  return(sum(x > 0))
}
otu.pruned$asv_occur_count = apply(otu.pruned, 1, ASVoccur)
summary(otu.pruned$asv_occur_count) # Minimum is 1
otu.highfreq = subset(otu.pruned, otu.pruned$asv_occur_count > 2)
summary(otu.highfreq$asv_occur_count) # Now minimum is 3
otu.highfreq = otu.highfreq[, -c(widthotu)]

## 4. Denoising
otu.clean <- mutate_all(otu.highfreq, funs(ifelse(. < 5, 0, .)))
ps_filt = phyloseq(sample_data(ps@sam_data),
                     tax_table(ps@tax_table),
                     otu_table(as.matrix(otu.clean), taxa_are_rows = TRUE))
ps_filt # 279 taxa, 132 samples
ps_filt@sam_data$sample_sums_filtered = sample_sums(ps_filt)
View(ps_filt@sam_data)

## Add % change in surface area growth measurements
growth_df_ali <- read.csv("growth_measurements_ali.csv")
growth_df_gurkiran <- read.csv("growth_measurements_gurkiran.csv")
growth_df_ali <- growth_df_ali %>%
  mutate(label = paste0(Origin, "_", Temperature..C.., "_", Well)) %>%
  select(label, Week.0.average.surface.area..in...)
growth_df_gurkiran <- growth_df_gurkiran %>%
  mutate(label = paste0(Origin, "_", Temperature..C.., "_", Well)) %>%
  select(label, Week.0.average.surface.area..in...)
growth_df <- bind_rows(growth_df_ali, growth_df_gurkiran)

meta <- data.frame(ps_filt@sam_data)
meta_merged <- left_join(meta, growth_df, by = "label") %>%
  mutate(percent_change_surfacearea = (Growth_surfacearea / Week.0.average.surface.area..in...)*100)
rownames(meta_merged) <- meta_merged$label

sample_data(ps_filt) <- sample_data(meta_merged)

ps_filt_gametophytes <- ps_filt
sample_data(ps_filt_gametophytes) <- data.frame(sample_data(ps_filt_gametophytes)) %>%
  filter(
    if_else(Experiment != "Gametophytebank", 
            !is.na(Treatment) & !is.na(Temperature), 
            TRUE)
  ) %>%
  mutate(Treatment = case_when(
    str_detect(Treatment, "antibiotic") ~ "antibiotic",
    .default = Treatment # Keeps original if it doesn't match
  )) %>%
  mutate(Temperature = as.character(Temperature)) %>% 
  mutate(Temperature = case_when(
    is.na(Temperature) ~ "",
    .default = Temperature # Keeps original if it doesn't match
  )) %>%
  mutate(Strain = case_when(
    Experiment == "Ali" ~ "EK Exp. 2",
    Experiment == "Gametophytebank" ~ "Gametophytebank",
    experiment == "Gam. Exp. 7 -16C" ~ "EK Exp. 1",
    experiment == "Gam. Exp. 7 -10C" ~ "EK Exp. 1",
    .default = Strain  # Keeps the original if it doesn't match
  )) %>%
  mutate(overall_treat = paste(Species, Treatment, Temperature, Strain, sep = " "))


ps_filt_nereocystis_gurkiran <- subset_samples(ps_filt, Experiment == "Gurkiran")
sample_data(ps_filt_nereocystis_gurkiran) <- data.frame(sample_data(ps_filt_nereocystis_gurkiran)) %>%
  filter(!is.na(Treatment), !is.na(Experiment)) %>%
  mutate(Treatment = str_replace(Treatment, "antibiotic 1", "antibiotic")) %>%
  mutate(overall_treat = paste(Treatment, Temperature, Strain, sep = " "))

saveRDS(ps_filt, "filtered_phyloseq.RDS")
saveRDS(ps_filt_gametophytes, "filtered_phyloseq_gametophytes.RDS")
saveRDS(ps_filt_nereocystis_gurkiran, "filtered_phyloseq_nereocystis_gurkiran.RDS")

################################################################## 
############################ RAREFYING ########################### 
################################################################## 

## Sample ID should be row names and ASV ID should be column names
View(otu.clean) # No
rarecurve(
  as.data.frame( 
    t( 
      as.matrix( 
        otu.clean))),
  step=50, cex=0.5, label=FALSE)
View(ps_filt@sam_data)
set.seed(5)
which(sample_sums(ps_filt_gametophytes) < 6000) # Lose 2 samples: 7_16_1_A6, EK_11_E2
which(sample_sums(ps_filt_gametophytes) < 20000) # Lose 31 samples
rarekelp_6000 <- rarefy_even_depth(ps_filt_gametophytes, sample.size = 6000) # Lose 2 samples
rarekelp_20000 <- rarefy_even_depth(ps_filt_gametophytes, sample.size = 20000) # Lose 31 samples
rarekelp_6000@sam_data$rare_sample_sums = sample_sums(rarekelp_6000)
summary(rarekelp_6000@sam_data$rare_sample_sums) # All are 6000
rarekelp_20000@sam_data$rare_sample_sums = sample_sums(rarekelp_20000)
summary(rarekelp_20000@sam_data$rare_sample_sums) # All are 20000

which(sample_sums(ps_filt_nereocystis_gurkiran) < 8000) # Lose 0 samples
which(sample_sums(ps_filt_nereocystis_gurkiran) < 20000) # Lose 6 samples
rarekelp_nereocystis_8000 <- rarefy_even_depth(ps_filt_nereocystis_gurkiran, sample.size = 8000) # Lose 40 taxa
rarekelp_nereocystis_20000 <- rarefy_even_depth(ps_filt_nereocystis_gurkiran, sample.size = 20000) # Lose 6 samples, 40 taxa
rarekelp_nereocystis_8000@sam_data$rare_sample_sums = sample_sums(rarekelp_nereocystis_8000)
summary(rarekelp_nereocystis_8000@sam_data$rare_sample_sums) # All are 8000
rarekelp_nereocystis_20000@sam_data$rare_sample_sums = sample_sums(rarekelp_nereocystis_20000)
summary(rarekelp_nereocystis_20000@sam_data$rare_sample_sums) # All are 20000

#saveRDS(rarekelp, "rarefied_phyloseq.RDS")
saveRDS(rarekelp_6000, "rarefied_phyloseq_gametophytes_6000reads.RDS")
saveRDS(rarekelp_20000, "rarefied_phyloseq_gametophytes_20000reads.RDS")
saveRDS(rarekelp_nereocystis_8000, "rarefied_phyloseq_nereocystis_gurkiran_8000reads.RDS")
saveRDS(rarekelp_nereocystis_20000, "rarefied_phyloseq_nereocystis_gurkiran_20000reads.RDS")
