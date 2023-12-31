####AIM 2B Beta Diversity Analysis####

###PART 1###

##Load Necessary Packages##
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

##setting a seed##
set.seed(1) 

##Load files##
load(file="Lab_Notebook/Phyloseq/dorms_rare_sheetwashfreq.RData")

##Adding combined sex and sheet washing frequency column##
sample_data(dorms_rare)$sex_sheetwashfreq <- paste(sample_data(dorms_rare)$sex, sample_data(dorms_rare)$sheetwashfreq_binned)

##Removing rows with "na" for sheetwashing frequency##
dorms_metadata_no_na <- sample_data(dorms_rare)

##Generating diversity metrics##
alphadiv <- estimate_richness(dorms_rare)
samp_dat_wdiv <- data.frame(dorms_metadata_no_na, alphadiv)

##Unweighted Unifrac##
unifrac_dm <- distance(dorms_rare, method="unifrac")
pcoa_unifrac <- ordinate(dorms_rare, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(dorms_rare, pcoa_unifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(dorms_rare, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT1/2B_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=8)



##Weighted Unifrac##
wu_dm <- distance(dorms_rare, method="wunifrac")
pcoa_wunifrac <- ordinate(dorms_rare, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(dorms_rare, pcoa_wunifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(dorms_rare, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT1/2B_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=8)



##Jaccard##
j_dm <- distance(dorms_rare, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(dorms_rare, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(dorms_rare, pcoa_jaccard, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(dorms_rare)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwashfreq, data=samp_dat_wdiv)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT1/2B_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=8)



##Bray Curtis##
bc_dm <- distance(dorms_rare, method="bray")
pcoa_bc <- ordinate(dorms_rare, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(dorms_rare, pcoa_bc, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(dorms_rare)), method="bray")
adonis2(dm_bray ~ sex_sheetwashfreq, data=samp_dat_wdiv)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT1/2B_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=8)
