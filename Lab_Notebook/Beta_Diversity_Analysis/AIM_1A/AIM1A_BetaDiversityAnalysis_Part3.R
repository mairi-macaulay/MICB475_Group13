####AIM 1A Beta Diversity Analysis PART 3####

##Load Necessary Packages##
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

##setting a seed##
set.seed(1) 

##Load files##
load(file="Lab_Notebook/Phyloseq/dorms_rare_sheetwashfreq.RData")

##Removing rows with "na" for sheetwashing frequency##
dorms_metadata_no_na <- sample_data(dorms_rare)

##Generating diversity metrics##
alphadiv <- estimate_richness(dorms_rare)
samp_dat_wdiv <- data.frame(dorms_metadata_no_na, alphadiv)

##Create new phyloseq with a metadata column that combines low and medium##
dorms_metadata_no_na$sheetwash_binned_combined_lowmed<- ifelse(dorms_metadata_no_na$sheetwashfreq_binned == "high", "high", "low/medium")
sample_data(dorms_rare) <- dorms_metadata_no_na

##Create new metadata with combined low and medium##
samp_dat_wdiv$sheetwash_binned_combined_lowmed <- ifelse(samp_dat_wdiv$sheetwashfreq_binned == "high", "high", "low/medium")

##Unweighted Unifrac##
unifrac_dm <- distance(dorms_rare, method="unifrac")
pcoa_unifrac <- ordinate(dorms_rare, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(dorms_rare, pcoa_unifrac, color = "sheetwash_binned_combined_lowmed") +
  labs(col="Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(dorms_rare, weighted=FALSE)
adonis2(unifrac_dm ~ sheetwash_binned_combined_lowmed, data=samp_dat_wdiv)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_3_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(dorms_rare, method="wunifrac")
pcoa_wunifrac <- ordinate(dorms_rare, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(dorms_rare, pcoa_unifrac, color = "sheetwash_binned_combined_lowmed") +
  labs(col="Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_unifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(dorms_rare, weighted=TRUE)
adonis2(wunifrac_dm ~ sheetwash_binned_combined_lowmed, data=samp_dat_wdiv)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_3_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(dorms_rare, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(dorms_rare, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(dorms_rare, pcoa_jaccard, color = "sheetwash_binned_combined_lowmed") + labs(col = "Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(dorms_rare)), method="jaccard")
adonis2(dm_jaccard ~ sheetwash_binned_combined_lowmed, data=samp_dat_wdiv)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_3_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(dorms_rare, method="bray")
pcoa_bc <- ordinate(dorms_rare, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(dorms_rare, pcoa_bc, color = "sheetwash_binned_combined_lowmed") + labs(col = "Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(dorms_rare)), method="bray")
adonis2(dm_bray ~ sheetwash_binned_combined_lowmed, data=samp_dat_wdiv)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_3_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)