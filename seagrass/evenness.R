library(phyloseq)
library(vegan)
library(ggplot2)
library(car)

## set wd

ps <- readRDS("../../Data/filtered_seagrass_mesocosm.rds")
ps <- subset_samples(ps, sample_type == "rhizome" & experiment_day == 20)

metadata <- as.data.frame(sample_data(ps))
otu <- t(as.data.frame(otu_table(ps)))

## Shannon index
metadata$shan = diversity(otu, index="shannon")
## Peilous evenness 
metadata$pielou = metadata$shan/(log(metadata$shan))

p <- ggplot(metadata, aes(x = treat, y = pielou, fill = treat)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_classic()
ggasave("Results/evenness_boxplot.png", plot = p)
p

# ANOVA
metadata <- data.frame(metadata)
a1 = aov(pielou~treat, data=metadata)
## use the summary function to see the output of the anova
summary(a1)

# sample size is equal as we know, 10 per treatment

# parametric or not:
plot(a1, which = 2)

# variance
leveneTest(a1)

# post hoc tukey
TukeyHSD(a1)
