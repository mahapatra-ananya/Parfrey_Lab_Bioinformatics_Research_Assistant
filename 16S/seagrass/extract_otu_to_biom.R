library(phyloseq)
library(biomformat);packageVersion("biomformat")

seagrass <- readRDS("../../Data/filtered_seagrass_mesocosm.RDS")

#otu <- t(as(otu_table(seagrass),"matrix")) # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE:
otu<-as(otu_table(seagrass),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"../../Data/otu_biom.biom")


######## running on select ASVs ######## 
selected_asvs <- c(
  "ASV7", "ASV12", "ASV22", "ASV26",
  "ASV30", "ASV33", "ASV55", "ASV68",
  "ASV130", "ASV163", "ASV169", "ASV197"
)
seagrass.subset <- prune_taxa(taxa_names(seagrass) %in% selected_asvs, seagrass)
otu.subset<-as(otu_table(seagrass.subset),"matrix")
otu_biom.subset<-make_biom(data=otu.subset)
write_biom(otu_biom.subset,"../../Data/otu_biom_selected_taxa.biom")