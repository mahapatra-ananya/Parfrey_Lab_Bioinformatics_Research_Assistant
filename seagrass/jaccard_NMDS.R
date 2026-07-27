# Load libraries
library(tidyverse)
library(phyloseq)
library(vegan)

# Set working directory

# Read in rarefied phyloseq object
seagrass <- readRDS("../../Data/rarefied_seagrass_mesocosm.RDS")

# Subset to only include measurements from the last day, from the rhizome samples
seagrass.subset <- subset_samples(seagrass, sample_type == "rhizome" & experiment_day == 20)

# Ordinate by Jaccard
seagrass.ord = ordinate(seagrass.subset, "NMDS","jaccard", binary = TRUE)

# Warning messages:
# 1: In decostand(x, "max", 2, na.rm = na.rm) :
# result contains NaN, perhaps due to impossible mathematical
# operation
#
# 2: In decostand(x, "total", 1, na.rm = na.rm) :
# result contains NaN, perhaps due to impossible mathematical
# operation
#
# 3: In metaMDS(veganifyOTU(physeq), distance, ...) :
# WA scores were not calculated due to missing values

# Check if any rows or columns are entirely 0
otutab <- as.data.frame((as.matrix(otu_table(seagrass.subset@otu_table))))
otutab$asv_abundance = rowSums(otutab)
min(otutab$asv_abundance)

# Remove ASVs that appear less than 0 times in the dataset
otu.pruned = subset(otutab, otutab$asv_abundance >= 1)

# Put back into phyloseq object
newseagrass = phyloseq(sample_data(seagrass.subset@sam_data),
                       tax_table(seagrass.subset@tax_table),
                       otu_table(as.matrix(otu.pruned), taxa_are_rows = TRUE))

# Ordinate by Jaccard
newseagrass.ord = ordinate(newseagrass, "NMDS","jaccard", binary = TRUE)

# Check fit
stressplot(newseagrass.ord)

# Check stress
newseagrass.ord
# Stress is 0.09 -> pretty good!

# Make NMDS plot
plot_ordination(newseagrass, newseagrass.ord, color = "treat") +
  stat_ellipse(aes(fill = treat), geom = "polygon", alpha = 0.1, linewidth = 0.5, show.legend = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(
    name = "Treatment", 
    values = c(
      "control_intact" = "#E15759", 
      "control_removed" = "#A0CBE8", 
      "sulfidic_intact" = "#499894", 
      "sulfidic_removed" = "#FABFD2"
    ), 
    labels = c(
      "control_intact" = "Control Intact", 
      "control_removed" = "Control Removed", 
      "sulfidic_intact" = "Sulfidic Intact", 
      "sulfidic_removed" = "Sulfidic Removed"
    ),
    aesthetics = c("color", "fill")
  ) + 
  theme_classic() + 
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )

# Warning message:
# In MASS::cov.trob(data[, vars], wt = weight * nrow(data)) :
# Probable convergence failure
# This means ellipses not fully reliable

# PERMANOVA:

# Calculate the distance within the data using Jaccard
distanceinfo <- phyloseq::distance(newseagrass, method = "jaccard")
# Extract metadata
sample_df <- data.frame(sample_data(newseagrass))
# Calculate the betadispersion within each treatment
region.dispers <- betadisper(distanceinfo, sample_df$treat)
## Betadispersion test to see if all treatments have the same betadispersion
beta.region = permutest(region.dispers)
beta.region

# Betadispersion test IS significant (p-value 0.041) so variance unequal
# among treatment groups (can see visually as well)
# Will still proceed with PERMANOVA due to lack of alternates currently

# Extract data frames from phyloseq object
otu.df = as.data.frame(t(as.matrix(newseagrass@otu_table)))
meta.df = as.data.frame(as.matrix(newseagrass@sam_data))

# Merge the data by Sample IDs into one data frame, using the sample ID
motu = merge(meta.df, otu.df, by = 0)

# Count column numbers
metacols = ncol(meta.df) + 1

# Run permanova test differences in community composition between regions
perm = adonis2(motu[, -c(1:metacols)] ~ treat, data = motu, method = "jaccard")
perm

# Results are significant (p-value 0.001)

# Post-hoc
pairwise = pairwise.adonis(motu[, -c(1:metacols)], motu$treat)
pairwise

# Same results as for bray-curtis: significant between all pairs except 
# SI vs SR, and CI vs CR
