####AIM 2B Beta Diversity Analysis####

###PART 4###

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

#Create new phylogenetic object for only low/high#
phylo_nomed <- subset_samples(dorms_rare, sex_sheetwashfreq %in% c("female high", "female low", "male high", "male low"))

##Create new metadata for only low/high#
samp_dat_wdiv_nomed <- subset(samp_dat_wdiv, sex_sheetwashfreq %in% c("female high", "female low", "male high", "male low"))



##Unweighted Unifrac##
unifrac_dm <- distance(phylo_nomed, method="unifrac")
pcoa_unifrac <- ordinate(phylo_nomed, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_nomed, pcoa_unifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_nomed, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_nomed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=8)



##Weighted Unifrac##
wu_dm <- distance(phylo_nomed, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_nomed, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_nomed, pcoa_wunifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_nomed, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_nomed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=8)



##Jaccard##
j_dm <- distance(phylo_nomed, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_nomed, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_nomed, pcoa_jaccard, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_nomed)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwashfreq, data=samp_dat_wdiv_nomed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=8)



##Bray Curtis##
bc_dm <- distance(phylo_nomed, method="bray")
pcoa_bc <- ordinate(phylo_nomed, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_nomed, pcoa_bc, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_nomed)), method="bray")
adonis2(dm_bray ~ sex_sheetwashfreq, data=samp_dat_wdiv_nomed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=8)



##Unweighted Unifrac WITH NO ELLIPSE##
unifrac_dm <- distance(phylo_nomed, method="unifrac")
pcoa_unifrac <- ordinate(phylo_nomed, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_nomed, pcoa_unifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups")
gg_pcoa_unifrac

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_unweighted_pcoa_noE.png"
       , gg_pcoa_unifrac
       , height=4, width=8)


##Weighted Unifrac WITH NO ELLIPSE##
wu_dm <- distance(phylo_nomed, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_nomed, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_nomed, pcoa_wunifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups")
gg_pcoa_wunifrac

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_weighted_pcoa_noE.png"
       , gg_pcoa_wunifrac
       , height=4, width=8)



##Jaccard WITH NO ELLIPSE##
j_dm <- distance(phylo_nomed, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_nomed, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_nomed, pcoa_jaccard, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups")
gg_pcoa_jaccard

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_jaccard_pcoa_noE.png"
       , gg_pcoa_jaccard
       , height=4, width=8)



##Bray Curtis WITH NO ELLIPSE BOLD##
bc_dm <- distance(phylo_nomed, method="bray")
pcoa_bc <- ordinate(phylo_nomed, method="PCoA", distance=bc_dm)
gg_pcoa_bc_noE <- plot_ordination(phylo_nomed, pcoa_bc, color = "sex_sheetwashfreq") + 
  labs(col = "Sheet Wash Frequency Sex Groups") +
  theme(text=element_text(size=16, face="bold"))
gg_pcoa_bc_noE

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_bc_pcoa_noE.png"
       , gg_pcoa_bc_noE
       , height=4, width=8)



##Bray Curtis WITH NO ELLIPSE##
bc_dm <- distance(phylo_nomed, method="bray")
pcoa_bc <- ordinate(phylo_nomed, method="PCoA", distance=bc_dm)
gg_pcoa_bc_noE <- plot_ordination(phylo_nomed, pcoa_bc, color = "sex_sheetwashfreq") + 
  labs(col = "Sex-Specific Sheet Washing Frequency Groups")
gg_pcoa_bc_noE

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT4/2B_4_bc_pcoa_noE_noB.png"
       , gg_pcoa_bc_noE
       , height=4, width=8)

