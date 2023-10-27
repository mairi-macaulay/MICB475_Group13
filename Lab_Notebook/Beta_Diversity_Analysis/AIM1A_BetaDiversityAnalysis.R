####AIM 1A Beta Diversity Analysis####

##Load Necessary Packages##
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)


##Load files##
dorms_metadata <- read_delim(file = "Lab_Notebook/metadata/dorms_metadata_updated.txt", delim = "\t")
load(file="Lab_Notebook/Phyloseq/dorms_rare_sheetwashfreq.RData")

##Removing rows with "na" for sheetwashing frequency##
dorms_metadata_no_na <- sample_data(dorms_rare)

##Generating diversity metrics##
alphadiv <- estimate_richness(dorms_rare)
samp_dat_wdiv <- data.frame(dorms_metadata_no_na, alphadiv)

##Unweighted Unifrac##
unifrac_dm <- distance(dorms_rare, method="unifrac")

pcoa_unifrac <- ordinate(dorms_rare, method="PCoA", distance=unifrac_dm)

gg_pcoa_unifrac <- plot_ordination(dorms_rare, pcoa_unifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency")
gg_pcoa_unifrac

#Statistical Analysis#
set.seed(1) #setting a seed like this ensures the stats test will be repeatable
unifrac_dm <- UniFrac(dorms_rare, weighted=FALSE)
adonis2(unifrac_dm ~ sheetwashfreq_binned, data=samp_dat_wdiv)


##Weighted Unifrac##
wu_dm <- distance(dorms_rare, method="wunifrac")

pcoa_wunifrac <- ordinate(dorms_rare, method="PCoA", distance=wu_dm)

gg_pcoa_wunifrac <- plot_ordination(dorms_rare, pcoa_wunifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency")
gg_pcoa_wunifrac

#Statistical Analysis#
wunifrac_dm <- UniFrac(dorms_rare, weighted=TRUE)
adonis2(wunifrac_dm ~ sheetwashfreq_binned, data=samp_dat_wdiv)











