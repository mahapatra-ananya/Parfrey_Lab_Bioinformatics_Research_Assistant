library(phyloseq)
library(dplyr)
library(vegan)
library(patchwork)
library(ggplot2)
library(car)
library(pairwiseAdonis)
library(ggplotify)
library(rstatix)

## set wd

kelp_ps <- readRDS("rarefied_phyloseq.RDS")
nereocystis_ps <- subset_samples(kelp_ps, Experiment == "Gurkiran")

gametophytes_ps_6000 <- readRDS("rarefied_phyloseq_gametophytes_6000reads.RDS")
gametophytes_ps_20000 <- readRDS("rarefied_phyloseq_gametophytes_20000reads.RDS")
nereocystis_ps_8000 <- readRDS("rarefied_phyloseq_nereocystis_gurkiran_8000reads.RDS")
nereocystis_ps_20000 <- readRDS("rarefied_phyloseq_nereocystis_gurkiran_20000reads.RDS")

gametophyte_ps_nonrarefied <- readRDS("filtered_phyloseq_gametophytes.RDS")
saccharina_ps <- subset_samples(gametophyte_ps_nonrarefied, Species == "Saccharina" & Experiment != "Gametophytebank")

gametophyte_colors <- c(
  # Nereocystis
  "Nereocystis control 16 CR"       = "red3",  
  "Nereocystis antibiotic 16 CR"    = "pink",  
  "Nereocystis control 10 CR"       = "blue4",  
  "Nereocystis antibiotic 10 CR"    = "lightblue", 
  "Nereocystis control 16 WS"       = "orange4",  
  "Nereocystis antibiotic 16 WS"    = "orange",  
  "Nereocystis control 10 WS"       = "green4",  
  "Nereocystis antibiotic 10 WS"    = "green",  
  # Saccharina
  "Saccharina control 16 EK Exp. 1"        = "#CC5500",  # purple
  "Saccharina antibiotic 16 EK Exp. 1"     = "#FFB347",  # orchid
  "Saccharina control 10 EK Exp. 1"        = "#008B8B",  # dark cyan
  "Saccharina antibiotic 10 EK Exp. 1"     = "#40E0D0",  # turquoise
  "Saccharina control 10 EK Exp. 2"       = "#4B0082",  # indigo
  "Saccharina antibiotic 10 EK Exp. 2"    = "#9370DB",  # medium purple
  "Saccharina   Gametophytebank"      = "#800080",  # brown
  "water   Gametophytebank"           = "#DA70D6"   # tan
)

nereocystis_colors <- c(
  "control 16 CR"    = "red3", # Dark Red
  "antibiotic 16 CR" = "pink", # Light Red
  "control 10 CR"    = "blue4", # Dark Blue
  "antibiotic 10 CR" = "lightblue", # Light Blue
  "control 16 WS"    = "orange4", # Dark Orange
  "antibiotic 16 WS" = "orange", # Light Orange
  "control 10 WS"    = "green4", # Dark Green
  "antibiotic 10 WS" = "green"  # Light Green
)

saccharina_colors <- c(
  "Saccharina control 16 EK Exp. 1"    = "red3", # Dark Red
  "Saccharina antibiotic 16 EK Exp. 1" = "pink", # Light Red
  "Saccharina control 10 EK Exp. 1"    = "blue4", # Dark Blue
  "Saccharina antibiotic 10 EK Exp. 1" = "lightblue", # Light Blue
  "Saccharina control 10 EK Exp. 2"    = "green4", # Dark Green
  "Saccharina antibiotic 10 EK Exp. 2" = "green"  # Light Green
)

####################################################################
######################### ALPHA DIVERSITY ##########################
####################################################################

calculate_alpha_div <- function(ps) {
  meta = as.data.frame(as.matrix(ps@sam_data))
  otu = t(as.matrix(ps@otu_table))
  meta$rich = specnumber(otu)
  meta$shan = diversity(otu, index="shannon")
  meta$pielou = meta$shan/(log(meta$rich))
  sample_data(ps) <- meta
  return(meta)
}

create_alpha_boxplots <- function(meta, treatment_var, treatment_str, alpha_metric, alpha_str, colors) {
  p <- ggplot(meta, aes(x=.data[[treatment_var]], y = .data[[alpha_metric]], fill=.data[[treatment_var]]))+
    geom_boxplot() +
    geom_jitter(width=0.1, height=0.1) +
    theme_bw() +
    labs(title = paste0(alpha_str, " Across ", treatment_str), 
         x=treatment_var, y=alpha_str) +
    theme(plot.title = element_text(size = 12), 
          legend.title = element_text(size = 10), 
          legend.text = element_text(size=8),
          axis.title = element_text(size=10),
          axis.text = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = colors)
}

create_alpha_plots <- function(meta, treatment_var, treatment_str, colors) {
  p1 <- create_alpha_boxplots(meta, treatment_var, treatment_str, "rich", "Richness", colors)
  p2 <- create_alpha_boxplots(meta, treatment_var, treatment_str, "pielou", "Evenness", colors)
  p3 <- create_alpha_boxplots(meta, treatment_var, treatment_str, "shan", "Shannon's Diversity Index", colors)
  p4 <- ggplot(meta, aes(x=rich, y = shan, color=.data[[treatment_var]]))+
    geom_point(alpha=0.7)+
    geom_smooth(method=lm, se=FALSE)+
    labs(title = paste0("Shannon Diversity vs. Richness at Different ", treatment_str), 
         x="Richness", y="Shannon's Diversity Index", color=treatment_var) +
    theme_bw() +
    theme(plot.title = element_text(size = 12), 
          legend.title = element_text(size = 10), 
          legend.text = element_text(size=8),
          axis.title = element_text(size=10),
          axis.text = element_text(size = 8)) +
    scale_color_manual(values = colors)
  p5 <- (p1+p2+p3) + plot_layout(guides = "collect")
  return(list("rich_box" = p1, "even_box" = p2, "shan_box" = p3, "rich_shan_dot" = p4, "combined" = p5))
}

check_variance <- function(meta, treatment_var) {
  print("Testing variance in richness: ")
  print(leveneTest(as.formula(paste0("rich ~", treatment_var)), data = meta))
  print("Testing variance in evenness: ")
  print(leveneTest(as.formula(paste0("pielou ~", treatment_var)), data = meta))
  print("Testing variance in shannon's diversity: ")
  print(leveneTest(as.formula(paste0("shan ~", treatment_var)), data = meta))
}

check_normality <- function(meta, alpha_metric, treatment_var) {
  anova_model <- aov(as.formula(paste0(alpha_metric, "~", treatment_var)), data = meta)
  plot(anova_model, which = 2)
  print(shapiro.test(residuals(anova_model)))
}

check_histogram <- function(meta, treatment_var, alpha_metric, alpha_str) {
  meta[[treatment_var]] <- factor(meta[[treatment_var]])
  treatment_levels <- levels(meta[[treatment_var]])
  par(mfrow = c(1, 2))
  hist(meta[[alpha_metric]][meta[[treatment_var]] == treatment_levels[1]],
       main = treatment_levels[1],
       xlab = paste0(alpha_str))
  hist(meta[[alpha_metric]][meta[[treatment_var]] == treatment_levels[2]],
       main = treatment_levels[2],
       xlab = paste0(alpha_str))
  par(mfrow = c(1, 1))
}

calculate_alpha_stat_parametric <- function(meta, treatment_var, alpha_metric, alpha_str, var_bool) {
  print(paste0("Comparing ", alpha_str, " :"))
  if (var_bool) {
    a <- aov(as.formula(paste0(alpha_metric, "~", treatment_var)), data = meta)
    print("ANOVA:")
    print(summary(a))
    print("Post-hoc Tukey:")
    TukeyHSD(a)
  } else {
    a <- oneway.test(as.formula(paste0(alpha_metric, "~", treatment_var)), data = meta, var.equal = FALSE) # Welch's ANOVA
    print("Welch's ANOVA")
    print(a)
    print("Post-hoc Games-Howell")
    print(games_howell_test(as.formula(paste0(alpha_metric, "~", treatment_var)), data = meta), n = Inf) # Games-Howell post-hoc
  }
}

calculate_alpha_stat_nonparametric <- function(meta, treatment_var, alpha_metric, alpha_str) {
  print(paste0("Comparing ", alpha_str, " :"))
  print("Kruskal-wallis:")
  print(kruskal.test(as.formula(paste0(alpha_metric, "~", treatment_var)), data = meta))
  print("Post-hoc Pairwise Wilcoxon:")
  print(pairwise.wilcox.test(meta[[alpha_metric]], meta[[treatment_var]], p.adjust.method = "BH", exact=FALSE))
}


calculate_alpha_stats_parametric <- function(meta, treatment_var, var_r, var_e, var_s) {
  calculate_alpha_stat_parametric(meta, treatment_var, "rich", "Richness", var_r)
  calculate_alpha_stat_parametric(meta, treatment_var, "pielou", "Evenness", var_e)
  calculate_alpha_stat_parametric(meta, treatment_var, "shan", "Shannon's Diversity", var_s)
}

calculate_alpha_stats_nonparametric <- function(meta, treatment_var) {
  calculate_alpha_stat_nonparametric(meta, treatment_var, "rich", "Richness")
  calculate_alpha_stat_nonparametric(meta, treatment_var, "pielou", "Evenness")
  calculate_alpha_stat_nonparametric(meta, treatment_var, "shan", "Shannon's Diversity")
}

## GURKIRAN NEREOCYSTIS
# 8000 reads
meta_nereocystis_8000 <- calculate_alpha_div(nereocystis_ps_8000)
alpha_plots_nereocystis_8000 <- create_alpha_plots(meta_nereocystis_8000, "overall_treat", "Treatment", nereocystis_colors)
combined_alpha_plot_nereocystis_8000 <- alpha_plots_nereocystis_8000[["combined"]]
combined_alpha_plot_nereocystis_8000
ggsave("alpha_plots_nereocystis_gurkiran_8000.png", combined_alpha_plot_nereocystis_8000)
check_variance(meta_nereocystis_8000, "overall_treat") # Levene's insignificant for all
check_normality(meta_nereocystis_8000, "rich", "overall_treat") # Possibly nonparametric
# check_histogram(meta_nereocystis_8000, "overall_treat", "rich", "Richness")
check_normality(meta_nereocystis_8000, "pielou", "overall_treat") # Normal (?)
# check_histogram(meta_nereocystis_8000, "overall_treat", "pielou", "Evenness")
check_normality(meta_nereocystis_8000, "shan", "overall_treat") # Normal
calculate_alpha_stat_nonparametric(meta_nereocystis_8000, "overall_treat", "rich", "Richness") # significant kruskal-wallis
calculate_alpha_stat_parametric(meta_nereocystis_8000, "overall_treat", "pielou", "Evenness", TRUE) # significant anova
calculate_alpha_stat_parametric(meta_nereocystis_8000, "overall_treat", "shan", "Shannon's Diversity", TRUE) # significant anova

# 20,000 reads
meta_nereocystis_20000 <- calculate_alpha_div(nereocystis_ps_20000)
alpha_plots_nereocystis_20000 <- create_alpha_plots(meta_nereocystis_20000, "overall_treat", "Treatment", nereocystis_colors)
combined_alpha_plot_nereocystis_20000 <- alpha_plots_nereocystis_20000[["combined"]]
combined_alpha_plot_nereocystis_20000
ggsave("alpha_plots_nereocystis_gurkiran_20000.png", combined_alpha_plot_nereocystis_20000)
check_variance(meta_nereocystis_20000, "overall_treat") # Levene's significant for richness
check_normality(meta_nereocystis_20000, "rich", "overall_treat") # Normal
# check_histogram(meta_nereocystis_20000, "overall_treat", "rich", "Richness")
check_normality(meta_nereocystis_20000, "pielou", "overall_treat") # Normal
# check_histogram(meta_nereocystis_20000, "overall_treat", "pielou", "Evenness")
check_normality(meta_nereocystis_20000, "shan", "overall_treat") # Normal
calculate_alpha_stat_parametric(meta_nereocystis_20000, "overall_treat", "rich", "Richness", FALSE) # significant welch's anova
calculate_alpha_stat_parametric(meta_nereocystis_20000, "overall_treat", "pielou", "Evenness", TRUE) # significant anova
calculate_alpha_stat_parametric(meta_nereocystis_20000, "overall_treat", "shan", "Shannon's Diversity", TRUE) # significant anova

###################################################################
######################### BETA DIVERSITY ##########################
###################################################################

calculate_and_plot_beta_ord <- function(ps, treatment_var, subset_var, colors, ord_bray, ord_jaccard) {
  sample_data(ps)[[treatment_var]] <- factor(sample_data(ps)[[treatment_var]])
  print("Stress for bray-curtis: ")
  print(ord_bray)
  ord_bray_p_ellipse <- plot_ordination(ps, ord_bray, color=treatment_var)+
    geom_point(size=3) +
    stat_ellipse() + theme_bw()+scale_color_manual(values = colors) +
    labs(title=paste0("Bray-curtis NMDS for ", subset_var))
  bray_df <- plot_ordination(ps, ord_bray, justDF = TRUE)
  agg_formula_bray <- reformulate(treatment_var, response = "cbind(NMDS1, NMDS2)")
  centroids_bray <- aggregate(agg_formula_bray, data = bray_df, FUN = mean)
  bray_build <- ggplot_build(ord_bray_p_ellipse)
  x_lim_bray <- bray_build$layout$panel_params[[1]]$x.range
  y_lim_bray <- bray_build$layout$panel_params[[1]]$y.range
  x_breaks_bray <- bray_build$layout$panel_params[[1]]$x$get_breaks()
  y_breaks_bray <- bray_build$layout$panel_params[[1]]$y$get_breaks()
  ord_bray_p_centroid <- ggplot(centroids_bray, 
                                aes(x = NMDS1, y = NMDS2, color = .data[[treatment_var]])) +
    geom_point(shape = 18, size = 5, show.legend = FALSE) +
    scale_x_continuous(limits = x_lim_bray, breaks = x_breaks_bray) +
    scale_y_continuous(limits = y_lim_bray, breaks = y_breaks_bray) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(title = paste0("Bray-curtis NMDS for ", subset_var),
         caption = paste("Stress =", round(ord_bray$stress, 4)))
  print("Stress for Jaccard: ")
  print(ord_jaccard)
  ord_jaccard_p_ellipse <- plot_ordination(ps, ord_jaccard, color=treatment_var)+
    geom_point(size=3) +
    stat_ellipse() + theme_bw()+ scale_color_manual(values = colors) +
    labs(title=paste0("Jaccard NMDS for ", subset_var))
  jaccard_df <- plot_ordination(ps, ord_jaccard, justDF = TRUE)
  agg_formula_jaccard <- reformulate(treatment_var, response = "cbind(NMDS1, NMDS2)")
  centroids_jaccard <- aggregate(agg_formula_jaccard, data = jaccard_df, FUN = mean)
  jaccard_build <- ggplot_build(ord_jaccard_p_ellipse)
  x_lim_jaccard <- jaccard_build$layout$panel_params[[1]]$x.range
  y_lim_jaccard <- jaccard_build$layout$panel_params[[1]]$y.range
  x_breaks_jaccard <- jaccard_build$layout$panel_params[[1]]$x$get_breaks()
  y_breaks_jaccard <- jaccard_build$layout$panel_params[[1]]$y$get_breaks()
  ord_jaccard_p_centroid <- ggplot(centroids_jaccard, 
                                aes(x = NMDS1, y = NMDS2, color = .data[[treatment_var]])) +
    geom_point(shape = 18, size = 5, show.legend = FALSE) +
    scale_x_continuous(limits = x_lim_jaccard, breaks = x_breaks_jaccard) +
    scale_y_continuous(limits = y_lim_jaccard, breaks = y_breaks_jaccard) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(title = paste0("Jaccard NMDS for ", subset_var),
         caption = paste("Stress =", round(ord_jaccard$stress, 4)))
  combined_p <- (ord_bray_p_ellipse + ord_jaccard_p_ellipse) / (ord_bray_p_centroid + ord_jaccard_p_centroid) +
    plot_layout(guides = "collect")
  return(combined_p)
}

check_variance_beta <- function(ps, treatment_var) {
  sample_df <- data.frame(sample_data(ps))
  distanceinfo_bray <- phyloseq::distance(ps, method = "bray")
  dispers_bray <- betadisper(distanceinfo_bray, sample_df[[treatment_var]])
  beta.dispers_bray = permutest(dispers_bray)
  print("Bray-curtis betadispersion: ")
  print(beta.dispers_bray)
  distanceinfo_jaccard <- phyloseq::distance(ps, method = "jaccard")
  dispers_jaccard <- betadisper(distanceinfo_jaccard, sample_df[[treatment_var]])
  beta.dispers_jaccard = permutest(dispers_jaccard)
  print("Jaccard betadispersion: ")
  print(beta.dispers_jaccard)
}

# 1-way permanova
# calculate_beta_stats <- function(ps, treatment_var) {
#   otu.df = as.data.frame(t(as.matrix(ps@otu_table)))
#   meta.df = as.data.frame(as.matrix(ps@sam_data))
#   motu = merge(meta.df, otu.df, by = 0)
#   metacols = ncol(meta.df) + 1
#   perm_bray = adonis2(as.formula(paste0("motu[, -c(1:metacols)] ~", treatment_var)), data = motu, method = "bray")
#   print("Bray-curtis PERMANOVA")
#   print(perm_bray)
#   pairwise_bray = pairwise.adonis(motu[, -c(1:metacols)], motu[[treatment_var]], sim.method = "bray")
#   print("Pairwise post-hoc test for bray-curtis")
#   print(pairwise_bray)
#   perm_jaccard = adonis2(as.formula(paste0("motu[, -c(1:metacols)] ~", treatment_var)), data = motu, method = "jaccard")
#   print("Jaccard PERMANOVA")
#   print(perm_jaccard)
#   pairwise_jaccard = pairwise.adonis(motu[, -c(1:metacols)], motu[[treatment_var]], sim.method = "jaccard")
#   print("Pairwise post-hoc test for jaccard")
#   print(pairwise_jaccard)
# }

# 3-factor sequential permanova with no interactions
calculate_beta_stats <- function(ps, treatment_var) {
  otu.df = as.data.frame(t(as.matrix(ps@otu_table)))
  meta.df = as.data.frame(as.matrix(ps@sam_data))
  motu = merge(meta.df, otu.df, by = 0)
  metacols = ncol(meta.df) + 1
  perm_bray = adonis2(motu[, -c(1:metacols)] ~ Strain + Treatment + Temperature, data = motu, method = "bray", by = "terms")
  print("Bray-curtis PERMANOVA")
  print(perm_bray)
  perm_jaccard = adonis2(motu[, -c(1:metacols)] ~ Strain + Treatment + Temperature, data = motu, method = "jaccard", by = "terms", binary = TRUE)
  print("Jaccard PERMANOVA")
  print(perm_jaccard)
}



set.seed(123)
## FULL GAMETOPHYTE DATASET
# 6000 reads threshold
ord_bray_gametophyte_6000 <- ordinate(gametophytes_ps_6000, "NMDS", "bray")
stressplot(ord_bray_gametophyte_6000)
ord_jaccard_gametophyte_6000 <- ordinate(gametophytes_ps_6000, "NMDS", "jaccard", binary=TRUE)
stressplot(ord_jaccard_gametophyte_6000)
beta_plot_gametophyte_6000 <- calculate_and_plot_beta_ord(gametophytes_ps_6000, 
                                                          "overall_treat", 
                                                          "Gametophyte 6000 reads",
                                                          gametophyte_colors,
                                                          ord_bray_gametophyte_6000,
                                                          ord_jaccard_gametophyte_6000)
print(beta_plot_gametophyte_6000)
ggsave("beta_plot_gametophyte_6000reads.png", beta_plot_gametophyte_6000)
check_variance_beta(gametophytes_ps_6000, "overall_treat")
calculate_beta_stats(gametophytes_ps_6000, "overall_treat")

# 20,000 reads threshold
ord_bray_gametophyte_20000 <- ordinate(gametophytes_ps_20000, "NMDS", "bray")
stressplot(ord_bray_gametophyte_20000)
ord_jaccard_gametophyte_20000 <- ordinate(gametophytes_ps_20000, "NMDS", "jaccard", binary=TRUE)
stressplot(ord_jaccard_gametophyte_20000)
beta_plot_gametophyte_20000 <- calculate_and_plot_beta_ord(gametophytes_ps_20000, 
                                                          "overall_treat", 
                                                          "Gametophyte 20000 reads",
                                                          gametophyte_colors,
                                                          ord_bray_gametophyte_20000,
                                                          ord_jaccard_gametophyte_20000)
print(beta_plot_gametophyte_20000)
ggsave("beta_plot_gametophyte_20000reads.png", beta_plot_gametophyte_20000)
check_variance_beta(gametophytes_ps_20000, "overall_treat")
calculate_beta_stats(gametophytes_ps_20000, "overall_treat")

## GURKIRAN NEREOCYSTIS
# 8000 reads threshold
ord_bray_nereocystis_8000 <- ordinate(nereocystis_ps_8000, "NMDS", "bray")
stressplot(ord_bray_nereocystis_8000)
ord_jaccard_nereocystis_8000 <- ordinate(nereocystis_ps_8000, "NMDS", "jaccard", binary=TRUE)
stressplot(ord_jaccard_nereocystis_8000)
beta_plot_nereocystis_8000 <- calculate_and_plot_beta_ord(nereocystis_ps_8000, 
                                                          "overall_treat", 
                                                          "Nereocystis 8000 reads",
                                                          nereocystis_colors,
                                                          ord_bray_nereocystis_8000,
                                                          ord_jaccard_nereocystis_8000)
print(beta_plot_nereocystis_8000)
ggsave("beta_plot_nereocystis_gurkiran_8000reads.png", beta_plot_nereocystis_8000)
check_variance_beta(nereocystis_ps_8000, "overall_treat")
calculate_beta_stats(nereocystis_ps_8000, "overall_treat")

# 20,000 reads threshold
ord_bray_nereocystis_20000 <- ordinate(nereocystis_ps_20000, "NMDS", "bray")
stressplot(ord_bray_nereocystis_20000)
ord_jaccard_nereocystis_20000 <- ordinate(nereocystis_ps_20000, "NMDS", "jaccard", binary=TRUE)
stressplot(ord_jaccard_nereocystis_20000)
beta_plot_nereocystis_20000 <- calculate_and_plot_beta_ord(nereocystis_ps_20000, 
                                                          "overall_treat", 
                                                          "Nereocystis 20000 reads",
                                                          nereocystis_colors,
                                                          ord_bray_nereocystis_20000,
                                                          ord_jaccard_nereocystis_20000)
print(beta_plot_nereocystis_20000)
ggsave("beta_plot_nereocystis_gurkiran_20000reads.png", beta_plot_nereocystis_20000)
check_variance_beta(nereocystis_ps_20000, "overall_treat")
calculate_beta_stats(nereocystis_ps_20000, "overall_treat")

######################################################################
######################### GROWTH COMPARISON ##########################
######################################################################

make_growth_plot <- function(df, x, y, title, x_label, y_label, colors) {
  p <- ggplot(df, 
         aes(x=.data[[x]], 
             y =.data[[y]], 
             fill=.data[[x]]))+
    geom_boxplot() +
    geom_jitter(width=0.1, height=0.1) +
    theme_bw() +
    labs(title = title, 
         x=x_label, y=y_label) +
    theme(plot.title = element_text(size = 12), 
          legend.title = element_text(size = 10), 
          legend.text = element_text(size=8),
          axis.title = element_text(size=10),
          axis.text = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = colors)
  return(p)
}

## Nereocystis (can use 8000 reads rarefied because no loss of samples)
meta_nereocystis_8000$percent_change_surfacearea <- as.numeric(meta_nereocystis_8000$percent_change_surfacearea)
meta_nereocystis_8000$Number_sporophytes_total <- as.numeric(meta_nereocystis_8000$Number_sporophytes_total)
meta_nereocystis_8000$Size_sporophytes_average_mm <- as.numeric(meta_nereocystis_8000$Size_sporophytes_average_mm)

## Nereocystis surface area growth
growth_plot_nereocystis_8000 <- make_growth_plot(meta_nereocystis_8000, 
                                                 "overall_treat", "percent_change_surfacearea",
                                                 "% Change in Growth (surface area inches) across treatments",
                                                 "Treatment", "% Change in Growth",
                                                 nereocystis_colors)
growth_plot_nereocystis_8000
ggsave("Growth_across_treatments_nereocystis_gurkiran.png", growth_plot_nereocystis_8000)
leveneTest(percent_change_surfacearea~overall_treat, data = meta_nereocystis_8000) # Insignificant
check_normality(meta_nereocystis_8000, "percent_change_surfacearea", "overall_treat") # Possibly nonparametric
calculate_alpha_stat_nonparametric(meta_nereocystis_8000, "overall_treat", "percent_change_surfacearea", "% Change in Growth") # significant kruskal-wallis
calculate_alpha_stat_parametric(meta_nereocystis_8000, "overall_treat", "percent_change_surfacearea", "% Change in Growth", TRUE) # significant anova

# Does increase in # sporophytes occur with an increase in % change surface area growth?
cor.test(meta_nereocystis_8000$percent_change_surfacearea, meta_nereocystis_8000$Number_sporophytes_total)
meta_nereocystis_control_10 <- meta_nereocystis_8000 %>% filter(overall_treat == "control 10 WS")
cor.test(meta_nereocystis_control_10$percent_change_surfacearea, meta_nereocystis_control_10$Number_sporophytes_total)
meta_nereocystis_control_16 <- meta_nereocystis_8000 %>% filter(overall_treat == "control 16 WS")
cor.test(meta_nereocystis_control_16$percent_change_surfacearea, meta_nereocystis_control_16$Number_sporophytes_total)
meta_nereocystis_numbersporophytes <- meta_nereocystis_8000 %>% filter(!is.na(Number_sporophytes_total))
ggplot(data = meta_nereocystis_numbersporophytes, aes(x = percent_change_surfacearea, y = Number_sporophytes_total, color = overall_treat)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(title = "Correlation between # sporophytes and % change in surface area growth",
       subtitle = "r = 0.15, p-val = 0.75")

## Saccharina
meta_saccharina <- sample_data(saccharina_ps)
meta_saccharina$percent_change_surfacearea <- as.numeric(meta_saccharina$percent_change_surfacearea)
meta_saccharina$Number_sporophytes_total <- as.numeric(meta_saccharina$Number_sporophytes_total)
meta_saccharina$Size_sporophytes_average_mm <- as.numeric(meta_saccharina$Size_sporophytes_average_mm)
meta_saccharina_surfacearea <- data.frame(meta_saccharina) %>% filter(!is.na(percent_change_surfacearea))
meta_saccharina_numbersporophytes <- data.frame(meta_saccharina) %>% filter(!is.na(Number_sporophytes_total))
meta_saccharina_sizesporophytes <- data.frame(meta_saccharina) %>% filter(!is.na(Size_sporophytes_average_mm))
desired_order_numbersporophytes <- c("Saccharina antibiotic 10 EK Exp. 2", "Saccharina control 10 EK Exp. 2", 
                                     "Saccharina antibiotic 10 EK Exp. 1", "Saccharina control 10 EK Exp. 1")
meta_saccharina_numbersporophytes$overall_treat <- 
  factor(meta_saccharina_numbersporophytes$overall_treat , levels = desired_order_numbersporophytes)

## Saccharina surface area growth
growth_plot_saccharina <- make_growth_plot(meta_saccharina_surfacearea,
                                           "overall_treat", "percent_change_surfacearea",
                                           "% Change in Growth (surface area inches) across treatments",
                                           "Treatment", "% Change in Growth",
                                           saccharina_colors)
growth_plot_saccharina
#ggsave("Growth_across_treatments_saccharina.png", growth_plot_saccharina)

## Saccharina number of sporophytes
numbersporophytes_plot_saccharina <- make_growth_plot(meta_saccharina_numbersporophytes,
                                                      "overall_treat", "Number_sporophytes_total",
                                                      "Number of sporophytes across treatments",
                                                      "Treatment", "# Sporophytes", saccharina_colors)




numbersporophytes_plot_saccharina 
#ggsave("Numbersporophytes_across_treatments_saccharina.png", numbersporophytes_plot_saccharina)

## Saccharina size of sporophytes
sizesporophytes_plot_saccharina <- make_growth_plot(meta_saccharina_sizesporophytes,
                                                    "overall_treat", "Size_sporophytes_average_mm",
                                                    "Size of sporophytes (mm) across treatments",
                                                    "Treatment", "Sporophyte size (mm)",
                                                    saccharina_colors)
sizesporophytes_plot_saccharina
#ggsave("Sizesporophytes_across_treatments_saccharina.png", sizesporophytes_plot_saccharina)

## Combined plot
combined_saccharina_growth_plot <- (growth_plot_saccharina + theme(legend.position = "none")) + 
  numbersporophytes_plot_saccharina + 
  (sizesporophytes_plot_saccharina + theme(legend.position = "none")) +
  plot_layout(guides = "collect")
combined_saccharina_growth_plot
ggsave("combined_growth_plot_saccharina.png", combined_saccharina_growth_plot)

leveneTest(percent_change_surfacearea~overall_treat, data = meta_saccharina_surfacearea) # Insignificant
check_normality(meta_saccharina_surfacearea, "percent_change_surfacearea", "overall_treat") # Nonparametric
wilcox.test(percent_change_surfacearea~overall_treat, data = meta_saccharina_surfacearea) # Significant

leveneTest(Number_sporophytes_total~overall_treat, data = meta_saccharina_numbersporophytes) # Insignificant
check_normality(meta_saccharina_numbersporophytes, "Number_sporophytes_total", "overall_treat") # Possibly nonparametric
calculate_alpha_stat_nonparametric(meta_saccharina_numbersporophytes, "overall_treat", "Number_sporophytes_total", "# Sporophytes") # Significant

leveneTest(Size_sporophytes_average_mm~overall_treat, data = meta_saccharina_sizesporophytes) # Insignificant
check_normality(meta_saccharina_sizesporophytes, "Size_sporophytes_average_mm", "overall_treat") # Nonparametric
wilcox.test(Size_sporophytes_average_mm~overall_treat, data = meta_saccharina_sizesporophytes) # Not significant

# Does increase in # sporophytes occur with an increase in % change surface area growth?
cor.test(meta_saccharina_surfacearea$percent_change_surfacearea, meta_saccharina_surfacearea$Number_sporophytes_total)
meta_saccharina_surfacearea_antibiotic <- meta_saccharina_surfacearea %>% filter(overall_treat == "Saccharina antibiotic 10 EK Exp. 2")
cor.test(meta_saccharina_surfacearea_antibiotic$percent_change_surfacearea, meta_saccharina_surfacearea_antibiotic$Number_sporophytes_total)
meta_saccharina_surfacearea_control <- meta_saccharina_surfacearea %>% filter(overall_treat == "Saccharina control 10 EK Exp. 2")
cor.test(meta_saccharina_surfacearea_control$percent_change_surfacearea, meta_saccharina_surfacearea_control$Number_sporophytes_total)
ggplot(data = meta_saccharina_surfacearea, aes(x = percent_change_surfacearea, y = Number_sporophytes_total, color = overall_treat)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(title = "Correlation between # sporophytes and % change in surface area growth",
       subtitle = "r = 0.07, p-val = 0.73")

## Check for all 12 samples
# Nereocystis
growth_df_gurkiran <- read.csv("growth_measurements_gurkiran.csv")
growth_df_gurkiran$Temperature..C.. <- as.factor(growth_df_gurkiran$Temperature..C..)
growth_df_gurkiran_mod <- growth_df_gurkiran %>%
  mutate(Treatment = case_when(
    str_detect(Treatment, "AB") ~ "antibiotic",
    str_detect(Treatment, "C") ~ "control",
    .default = Treatment # Keeps original if it doesn't match
  )) %>%
  mutate(Temperature..C.. = case_when(
    Temperature..C.. == "11" ~ "10",
    .default = Temperature..C.. # Keeps original if it doesn't match
  )) %>%
  mutate(overall_treat = paste(Treatment, Temperature..C.., Origin, sep = " ")) %>%
  mutate(label = paste0(Origin, "_", Temperature..C.., "_", Well)) %>%
  mutate(percent_change_surfacearea = (Change.in.surface.area.coverage / Week.0.average.surface.area..in...)*100) %>%
  select(label, Treatment, Temperature..C.., Origin, overall_treat, percent_change_surfacearea) %>%
  filter(label != "_NA_")

full_growth_plot_nereocystis <- make_growth_plot(growth_df_gurkiran_mod,
                                                 "overall_treat", "percent_change_surfacearea",
                                                 "% Change in Growth (surface area inches) across treatments",
                                                 "Treatment", "% Change in Growth",
                                                 nereocystis_colors)
full_growth_plot_nereocystis
ggsave("full_growth_plot_nereocystis.png", full_growth_plot_nereocystis)
leveneTest(percent_change_surfacearea~overall_treat, data = growth_df_gurkiran_mod) # Significant
check_normality(growth_df_gurkiran_mod, "percent_change_surfacearea", "overall_treat") # Possibly nonparametric... looks parametric mostly
#calculate_alpha_stat_nonparametric(growth_df_gurkiran_mod, "overall_treat", "percent_change_surfacearea", "% Change in Growth") # significant kruskal-wallis
calculate_alpha_stat_parametric(growth_df_gurkiran_mod, "overall_treat", "percent_change_surfacearea", "% Change in Growth", FALSE) # significant anova

# Saccharina
saccharina_colors_full <- c(
  "Saccharina Control 16 EK"    = "red3", # Dark Red
  "Saccharina antibiotic 16 EK" = "pink", # Light Red
  "Saccharina Control 10 EK"    = "blue4", # Dark Blue
  "Saccharina antibiotic 10 EK" = "lightblue", # Light Blue
  "Saccharina Control 16 AB"    = "orange4", # Dark Orange
  "Saccharina antibiotic 16 AB" = "orange", # Light Orange
  "Saccharina Control 10 AB"    = "green4", # Dark Green
  "Saccharina antibiotic 10 AB" = "green"  # Light Green
)
growth_df_ali <- read.csv("growth_measurements_ali.csv")
growth_df_ali$Temperature..C.. <- as.factor(growth_df_ali$Temperature..C..)
growth_df_ali_mod <- growth_df_ali %>%
  filter(Species != "Nereocystis") %>%
  mutate(Treatment = case_when(
    str_detect(Treatment, "Antibiotic") ~ "antibiotic",
    .default = Treatment # Keeps original if it doesn't match
  )) %>%
  mutate(Temperature..C.. = case_when(
    Temperature..C.. == "11" ~ "10",
    .default = Temperature..C.. # Keeps original if it doesn't match
  )) %>%
  mutate(overall_treat = paste(Species, Treatment, Temperature..C.., Origin, sep = " ")) %>%
  mutate(label = paste0(Origin, "_", Temperature..C.., "_", Well)) %>%
  mutate(percent_change_surfacearea = (Change.in.surface.area.coverage / Week.0.average.surface.area..in...)*100) %>%
  select(label, Treatment, Temperature..C.., Origin, overall_treat, percent_change_surfacearea)
full_growth_plot_saccharina <- make_growth_plot(growth_df_ali_mod,
                                                 "overall_treat", "percent_change_surfacearea",
                                                 "% Change in Growth (surface area inches) across treatments",
                                                 "Treatment", "% Change in Growth",
                                                 saccharina_colors_full)
full_growth_plot_saccharina
ggsave("full_growth_plot_saccharina.png", full_growth_plot_saccharina)

growth_df_ali_10C <- growth_df_ali_mod %>%
  filter(Temperature..C.. == "10")
leveneTest(percent_change_surfacearea~overall_treat, data = growth_df_ali_10C) # Insignificant
check_normality(growth_df_ali_10C, "percent_change_surfacearea", "overall_treat") # Nonparametric
calculate_alpha_stat_nonparametric(growth_df_ali_10C, "overall_treat", "percent_change_surfacearea", "% Change in Growth")





#################### old code for subsetting nereocystis further ####################

# ASVoccur = function(x) {
#   return(sum(x > 0))
# }

# clean_subset_ps <- function(ps) {
#   otutab <- as.data.frame(as.matrix(otu_table(ps@otu_table)))
#   print(head(otutab))
#   otutab$asv_abundance = rowSums(otutab)
#   print(min(otutab$asv_abundance))
#   otu.pruned = subset(otutab, otutab$asv_abundance >= 1)
#   print(min(otu.pruned$asv_abundance))
#   widthotu = ncol(otu.pruned)
#   otu.pruned = otu.pruned[, -c(widthotu)]
#   # otu.pruned$asv_occur_count = apply(otu.pruned, 1, ASVoccur)
#   # print(summary(otu.pruned$asv_occur_count))
#   # return(otu.pruned)
# }
# 
# nereocystis_control <- subset_samples(nereocystis_ps, Treatment == "control")
# otu_control <- clean_subset_ps(nereocystis_control)
# otu_table(nereocystis_control) <- otu_table(otu_control, taxa_are_rows = TRUE)
# nereocystis_antib <- subset_samples(nereocystis_ps, Treatment != "control")
# otu_antib <- clean_subset_ps(nereocystis_antib)
# otu_table(nereocystis_antib) <- otu_table(otu_antib, taxa_are_rows = TRUE)
# nereocystis_10 <- subset_samples(nereocystis_ps, Temperature == 10)
# otu_10 <- clean_subset_ps(nereocystis_10)
# otu_table(nereocystis_10) <- otu_table(otu_10, taxa_are_rows = TRUE)
# nereocystis_16 <- subset_samples(nereocystis_ps, Temperature == 16)
# otu_16 <- clean_subset_ps(nereocystis_16)
# otu_table(nereocystis_16) <- otu_table(otu_16, taxa_are_rows = TRUE)

# check_normality <- function(meta, alpha_metric, treatment_var) {
#   t_model <- lm(as.formula(paste0(alpha_metric, "~", treatment_var)), data = meta)
#   plot(t_model, which = 2)
#   print(shapiro.test(residuals(t_model)))
# }
# calculate_alpha_stat_parametric <- function(meta, treatment_var, alpha_metric, alpha_str, var_bool) {
#   print(paste0("Comparing ", alpha_str, " :"))
#   print(t.test(as.formula(paste0(alpha_metric, "~", treatment_var)), data = meta, var.equal = var_bool))
# }
# calculate_alpha_stat_nonparametric <- function(meta, treatment_var, alpha_metric, alpha_str) {
#   print(paste0("Comparing ", alpha_str, " :"))
#   print(wilcox.test(as.formula(paste0(alpha_metric, "~", treatment_var)), data = meta))
# }

# ## CONTROL MICROBIAL TREATMENT
# meta_control <- calculate_alpha_div(nereocystis_control)
# alpha_plots_control <- create_alpha_plots(meta_control, "Temperature", "Temperatures")
# combined_alpha_plot_control <- alpha_plots_control[["combined"]]
# combined_alpha_plot_control
# ggsave("alpha_plots_controlmicrobe.png", combined_alpha_plot_control)
# check_variance(meta_control, "Temperature") # Levene's insignificant for all 3 alpha metrics
# check_normality(meta_control, "rich", "Temperature")
# check_histogram(meta_control, "Temperature", "rich", "Richness")
# check_normality(meta_control, "pielou", "Temperature")
# check_histogram(meta_control, "Temperature", "pielou", "Evenness")
# check_normality(meta_control, "shan", "Temperature")
# check_histogram(meta_control, "Temperature", "shan", "Shannon's Diversity")
# calculate_alpha_stats_nonparametric(meta_control, "Temperature") # Will use nonparam for all
# #calculate_alpha_stats_parametric(meta_control, "Temperature", TRUE, TRUE, TRUE) # also not significant
# 
# ## ANTIBIOTICS MICROBIAL TREATMENT
# meta_antib <- calculate_alpha_div(nereocystis_antib)
# alpha_plots_antib <- create_alpha_plots(meta_antib, "Temperature", "Temperatures")
# combined_alpha_plot_antib <- alpha_plots_antib[["combined"]]
# combined_alpha_plot_antib
# ggsave("alpha_plots_antibioticsmicrobe.png", combined_alpha_plot_antib)
# check_variance(meta_antib, "Temperature") # Levene's insignificant for all 3 alpha metrics
# check_normality(meta_antib, "rich", "Temperature")
# check_histogram(meta_antib, "Temperature", "rich", "Richness")
# check_normality(meta_antib, "pielou", "Temperature")
# check_histogram(meta_antib, "Temperature", "pielou", "Evenness")
# check_normality(meta_antib, "shan", "Temperature")
# check_histogram(meta_antib, "Temperature", "shan", "Shannon's Diversity")
# calculate_alpha_stats_nonparametric(meta_antib, "Temperature") # Will use nonparam for all
# 
# ## 10C 
# meta_10 <- calculate_alpha_div(nereocystis_10)
# alpha_plots_10 <- create_alpha_plots(meta_10, "Treatment", "Microbial Treatment")
# combined_alpha_plot_10 <- alpha_plots_10[["combined"]]
# combined_alpha_plot_10
# ggsave("alpha_plots_10C.png", combined_alpha_plot_10)
# check_variance(meta_10, "Treatment") # Levene's significant for evenness
# check_normality(meta_10, "rich", "Treatment")
# check_histogram(meta_10, "Treatment", "rich", "Richness")
# check_normality(meta_10, "pielou", "Treatment")
# check_histogram(meta_10, "Treatment", "pielou", "Evenness")
# check_normality(meta_10, "shan", "Treatment")
# check_histogram(meta_10, "Treatment", "shan", "Shannon's Diversity")
# calculate_alpha_stat_nonparametric(meta_10, "Treatment", "rich", "Richness")
# calculate_alpha_stat_parametric(meta_10, "Treatment", "pielou", "Evenness", FALSE)
# calculate_alpha_stat_parametric(meta_10, "Treatment", "shan", "Shannon's Diversity", TRUE)
# 
# ## 16C
# meta_16 <- calculate_alpha_div(nereocystis_16)
# alpha_plots_16 <- create_alpha_plots(meta_16, "Treatment", "Microbial Treatment")
# combined_alpha_plot_16 <- alpha_plots_16[["combined"]]
# combined_alpha_plot_16
# ggsave("alpha_plots_16C.png", combined_alpha_plot_16)
# check_variance(meta_16, "Treatment") # Levene's insignificant for all 3 alpha metrics
# check_normality(meta_16, "rich", "Treatment")
# check_histogram(meta_16, "Treatment", "rich", "Richness")
# check_normality(meta_16, "pielou", "Treatment")
# check_histogram(meta_16, "Treatment", "pielou", "Evenness")
# check_normality(meta_16, "shan", "Treatment")
# check_histogram(meta_16, "Treatment", "shan", "Shannon's Diversity")
# calculate_alpha_stat_parametric(meta_16, "Treatment", "rich", "Richness", TRUE)
# calculate_alpha_stat_nonparametric(meta_16, "Treatment", "pielou", "Evenness")
# calculate_alpha_stat_nonparametric(meta_16, "Treatment", "shan", "Shannon's Diversity")

# ## CONTROL MICROBIAL TREATMENT
# beta_plot_control <- calculate_and_plot_beta_ord(nereocystis_control, "Temperature", "control microbe")
# beta_plot_control[["stress_bray"]]
# beta_plot_control[["stress_jaccard"]]
# beta_plot_control[["combined"]]
# ggsave("beta_plots_control.png", beta_plot_control[["combined"]])
# check_variance_beta(nereocystis_control, "Temperature") # Not significant for either
# calculate_beta_stats(nereocystis_control, "Temperature")
# 
# ## ANTIBIOTICS MICROBIAL TREATMENT
# beta_plot_antib <- calculate_and_plot_beta_ord(nereocystis_antib, "Temperature", "antibiotics")
# beta_plot_antib[["stress_bray"]]
# beta_plot_antib[["stress_jaccard"]]
# beta_plot_antib[["combined"]]
# ggsave("beta_plots_antib.png", beta_plot_antib[["combined"]])
# check_variance_beta(nereocystis_antib, "Temperature") # Not significant for either
# calculate_beta_stats(nereocystis_antib, "Temperature")
# 
# ## 10C 
# beta_plot_10 <- calculate_and_plot_beta_ord(nereocystis_10, "Treatment", "10C")
# beta_plot_10[["stress_bray"]]
# beta_plot_10[["stress_jaccard"]]
# beta_plot_10[["combined"]]
# ggsave("beta_plots_10.png", beta_plot_10[["combined"]])
# check_variance_beta(nereocystis_10, "Treatment") # Not significant for either
# calculate_beta_stats(nereocystis_10, "Treatment")
# 
# ## 16C
# beta_plot_16 <- calculate_and_plot_beta_ord(nereocystis_16, "Treatment", "16C")
# beta_plot_16[["stress_bray"]]
# beta_plot_16[["stress_jaccard"]]
# beta_plot_16[["combined"]] 
# ggsave("beta_plots_16.png", beta_plot_16[["combined"]])
# check_variance_beta(nereocystis_16, "Treatment") # Not significant for either
# calculate_beta_stats(nereocystis_10, "Treatment")
