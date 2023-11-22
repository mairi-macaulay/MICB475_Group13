####AIM 2B Beta Diversity Analysis####

###PART 3###

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

##Add a column that combines low and medium##
dorms_metadata_no_na$sheetwash_binned_combined_lowmed<- ifelse(dorms_metadata_no_na$sheetwashfreq_binned == "high", "high", "low/medium")
sample_data(dorms_rare) <- dorms_metadata_no_na

##Adding combined sex and sheet washing frequency column##
dorms_metadata_no_na$sex_sheetwash_lowmedcomb <- paste(sample_data(dorms_rare)$sex, sample_data(dorms_rare)$sheetwash_binned_combined_lowmed)
sample_data(dorms_rare) <- dorms_metadata_no_na

##Generating diversity metrics##
alphadiv <- estimate_richness(dorms_rare)
samp_dat_wdiv <- data.frame(dorms_metadata_no_na, alphadiv)

##Create new phyloseq for each condition##
#female low/medium vs female high#
phylo_f_lowmedvshigh <- subset_samples(dorms_rare, sex_sheetwash_lowmedcomb %in% c("female high", "female low/medium"))
#male low/medium vs male high#
phylo_m_lowmedvshigh <- subset_samples(dorms_rare, sex_sheetwash_lowmedcomb %in% c("male high", "male low/medium"))
#female low/medium vs male low/medium#
phylo_mlowmedvsflowmed <- subset_samples(dorms_rare, sex_sheetwash_lowmedcomb %in% c("male low/medium", "female low/medium"))
#female high vs male high#
phylo_mhighvsfhigh <- subset_samples(dorms_rare, sex_sheetwash_lowmedcomb %in% c("male high", "female high"))

##Create new metadata for each condition##
#female low/medium vs female high#
samp_dat_wdiv_f_lowmedvshigh <- subset(samp_dat_wdiv, sex_sheetwash_lowmedcomb %in% c("female high", "female low/medium"))
#male low/medium vs male high#
samp_dat_m_lowmedvshigh <- subset(samp_dat_wdiv, sex_sheetwash_lowmedcomb %in% c("male high", "male low/medium"))
#female low/medium vs male low/medium#
samp_dat_wdiv_mlowmedvsflowmed <- subset(samp_dat_wdiv, sex_sheetwash_lowmedcomb %in% c("male low/medium", "female low/medium"))
#female high vs male high#
samp_dat_wdiv_mhighvsfhigh <- subset(samp_dat_wdiv, sex_sheetwash_lowmedcomb %in% c("female high", "male high"))


###FEMALE LOW/MEDIUM vs FEMALE HIGH###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_f_lowmedvshigh, method="unifrac")
pcoa_unifrac <- ordinate(phylo_f_lowmedvshigh, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_f_lowmedvshigh, pcoa_unifrac, color = "sex_sheetwash_lowmedcomb") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_f_lowmedvshigh, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_f_lowmedvshigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_F_highvslowmed_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_f_lowmedvshigh, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_f_lowmedvshigh, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_f_lowmedvshigh, pcoa_wunifrac, color = "sex_sheetwash_lowmedcomb") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_f_lowmedvshigh, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_f_lowmedvshigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_F_highvslowmed_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_f_lowmedvshigh, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_f_lowmedvshigh, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_f_lowmedvshigh, pcoa_jaccard, color = "sex_sheetwash_lowmedcomb") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_f_lowmedvshigh)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_f_lowmedvshigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_F_highvslowmed_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_f_lowmedvshigh, method="bray")
pcoa_bc <- ordinate(phylo_f_lowmedvshigh, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_f_lowmedvshigh, pcoa_bc, color = "sex_sheetwash_lowmedcomb") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_f_lowmedvshigh)), method="bray")
adonis2(dm_bray ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_f_lowmedvshigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_F_highvslowmed_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)



###MALE LOW/MEDIUM vs MALE HIGH###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_m_lowmedvshigh, method="unifrac")
pcoa_unifrac <- ordinate(phylo_m_lowmedvshigh, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_m_lowmedvshigh, pcoa_unifrac, color = "sex_sheetwash_lowmedcomb") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_m_lowmedvshigh, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwash_lowmedcomb, data=samp_dat_m_lowmedvshigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_M_highvslowmed_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_m_lowmedvshigh, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_m_lowmedvshigh, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_m_lowmedvshigh, pcoa_wunifrac, color = "sex_sheetwash_lowmedcomb") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_m_lowmedvshigh, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwash_lowmedcomb, data=samp_dat_m_lowmedvshigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_M_highvslowmed_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_m_lowmedvshigh, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_m_lowmedvshigh, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_m_lowmedvshigh, pcoa_jaccard, color = "sex_sheetwash_lowmedcomb") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_m_lowmedvshigh)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwash_lowmedcomb, data=samp_dat_m_lowmedvshigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_M_highvslowmed_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_m_lowmedvshigh, method="bray")
pcoa_bc <- ordinate(phylo_m_lowmedvshigh, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_m_lowmedvshigh, pcoa_bc, color = "sex_sheetwash_lowmedcomb") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_m_lowmedvshigh)), method="bray")
adonis2(dm_bray ~ sex_sheetwash_lowmedcomb, data=samp_dat_m_lowmedvshigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_M_highvslowmed_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)



###MALE LOW/MEDIUM vs FEMALE LOW/MEDIUM###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_mlowmedvsflowmed, method="unifrac")
pcoa_unifrac <- ordinate(phylo_mlowmedvsflowmed, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_mlowmedvsflowmed, pcoa_unifrac, color = "sex_sheetwash_lowmedcomb") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_mlowmedvsflowmed, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_mlowmedvsflowmed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_lowmed_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_mlowmedvsflowmed, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_mlowmedvsflowmed, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_mlowmedvsflowmed, pcoa_wunifrac, color = "sex_sheetwash_lowmedcomb") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_mlowmedvsflowmed, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_mlowmedvsflowmed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_lowmed_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_mlowmedvsflowmed, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_mlowmedvsflowmed, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_mlowmedvsflowmed, pcoa_jaccard, color = "sex_sheetwash_lowmedcomb") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_mlowmedvsflowmed)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_mlowmedvsflowmed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_lowmed_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_mlowmedvsflowmed, method="bray")
pcoa_bc <- ordinate(phylo_mlowmedvsflowmed, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_mlowmedvsflowmed, pcoa_bc, color = "sex_sheetwash_lowmedcomb") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_mlowmedvsflowmed)), method="bray")
adonis2(dm_bray ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_mlowmedvsflowmed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_lowmed_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)




###MALE HIGH vs FEMALE HIGH###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_mhighvsfhigh, method="unifrac")
pcoa_unifrac <- ordinate(phylo_mhighvsfhigh, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_mhighvsfhigh, pcoa_unifrac, color = "sex_sheetwash_lowmedcomb") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_mhighvsfhigh, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_mhighvsfhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_high_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_mhighvsfhigh, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_mhighvsfhigh, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_mhighvsfhigh, pcoa_wunifrac, color = "sex_sheetwash_lowmedcomb") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_mhighvsfhigh, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_mhighvsfhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_high_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_mhighvsfhigh, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_mhighvsfhigh, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_mhighvsfhigh, pcoa_jaccard, color = "sex_sheetwash_lowmedcomb") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_mhighvsfhigh)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_mhighvsfhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_high_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_mhighvsfhigh, method="bray")
pcoa_bc <- ordinate(phylo_mhighvsfhigh, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_mhighvsfhigh, pcoa_bc, color = "sex_sheetwash_lowmedcomb") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_mhighvsfhigh)), method="bray")
adonis2(dm_bray ~ sex_sheetwash_lowmedcomb, data=samp_dat_wdiv_mhighvsfhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT3/2B_3_high_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)