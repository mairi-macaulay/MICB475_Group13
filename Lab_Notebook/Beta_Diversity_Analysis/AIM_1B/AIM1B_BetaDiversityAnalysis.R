####AIM 1B Beta Diversity Analysis####

##Load Necessary Packages##
library(phyloseq)
library(ape)
library(tidyverse)
library(picante) 

##Set a seed##
set.seed(2)

##Load files##
load(file="Lab_Notebook/Phyloseq/dorms_rare_showerrecency.RData")

##Removing rows with "na" for sheetwashing frequency##
dorms_metadata_no_na <- sample_data(dorms_rare)

##Generating diversity metrics##
alphadiv <- estimate_richness(dorms_rare)
samp_dat_wdiv <- data.frame(dorms_metadata_no_na, alphadiv)


##Unweighted Unifrac##
unifrac_dm <- distance(dorms_rare, method="unifrac")
pcoa_unifrac <- ordinate(dorms_rare, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(dorms_rare, pcoa_unifrac, color = "last_shower_binned") +
  labs(col="Shower Recency") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
unifrac_dm <- UniFrac(dorms_rare, weighted=FALSE)
adonis2(unifrac_dm ~ last_shower_binned, data=samp_dat_wdiv)

#Save plots#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1B/1B_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)



##Weighted Unifrac##
wu_dm <- distance(dorms_rare, method="wunifrac")
pcoa_wunifrac <- ordinate(dorms_rare, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(dorms_rare, pcoa_wunifrac, color = "last_shower_binned") +
  labs(col="Shower Recency") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
wunifrac_dm <- UniFrac(dorms_rare, weighted=TRUE)
adonis2(wunifrac_dm ~ last_shower_binned, data=samp_dat_wdiv)

#Save plots#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1B/1B_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)



##Jaccard##
j_dm <- distance(dorms_rare, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(dorms_rare, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(dorms_rare, pcoa_jaccard, color = "last_shower_binned") + labs(col = "Shower Recency") + stat_ellipse()
pcoa_jaccard

#Jaccard diversity PERMANOVA Test#
dm_jaccard <- vegdist(t(otu_table(dorms_rare)), method="jaccard")
adonis2(dm_jaccard ~ last_shower_binned, data=samp_dat_wdiv)

#Save plots#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1B/1B_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)



##Bray Curtis##
bc_dm <- distance(dorms_rare, method="bray")
pcoa_bc <- ordinate(dorms_rare, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(dorms_rare, pcoa_bc, color = "last_shower_binned") + labs(col = "Shower Recency") + stat_ellipse()

#Bray Curtic diversity PERMANOVA Test#
dm_bray <- vegdist(t(otu_table(dorms_rare)), method="bray")
adonis2(dm_bray ~ last_shower_binned, data=samp_dat_wdiv)

#Save plots#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1B/1B_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)


