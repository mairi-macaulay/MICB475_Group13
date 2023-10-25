####AIM 1A Beta Diversity Analysis####

##Load Necessary Packages##
library(phyloseq)
library(ape)
library(tidyverse)
library(picante) 

dorms_metadata <- read_delim(file = "dorms_metadata_updated.txt", delim = "\t")

##Filter out "NA" values##
dorms_rare <- subset_samples(dorms_rare, sheetwashfreq_binned != "NA")

##Unweighted Unifrac##
unifrac_dm <- distance(dorms_rare, method="unifrac")

pcoa_unifrac <- ordinate(dorms_rare, method="PCoA", distance=unifrac_dm)

plot_ordination(dorms_rare, pcoa_unifrac, color = "sheetwashfreq_binned")

gg_pcoa_unifrac <- plot_ordination(dorms_rare, pcoa_unifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency")
gg_pcoa_unifrac

unifrac_dm <- UniFrac(dorms_rare, weighted=FALSE)
adonis2(unifrac_dm ~ sheet.wash.frequency, data=)


##Weighted Unifrac##
wunifrac_dm <- distance(dorms_rare, method="wunifrac")

pcoa_wunifrac <- ordinate(dorms_rare, method="PCoA", distance=wunifrac_dm)

plot_ordination(dorms_rare, pcoa_wunifrac, color = "host_subject_id", shape="sheetwashfreq_binned")

gg_pcoa_wunifrac <- plot_ordination(dorms_rare, pcoa_wunifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency")
gg_pcoa_wunifrac

wunifrac_dm <- UniFrac(dorms_rare, weighted=TRUE)
adonis2(wunifrac_dm ~ sheet.wash.frequency, data=)










