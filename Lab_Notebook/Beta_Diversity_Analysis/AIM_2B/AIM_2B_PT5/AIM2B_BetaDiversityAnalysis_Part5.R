####AIM 2B Beta Diversity Analysis####

###PART 5###

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

##Create new phyloseq for each condition##
#female low vs female high#
phylo_f_lowvshigh <- subset_samples(dorms_rare, sex_sheetwashfreq %in% c("female high", "female low"))
#male low vs male high#
phylo_m_lowvshigh <- subset_samples(dorms_rare, sex_sheetwashfreq %in% c("male high", "male low"))
#female low vs male low#
phylo_mlowvsflow <- subset_samples(dorms_rare, sex_sheetwashfreq %in% c("male low", "female low"))
#female high vs male high#
phylo_mhighvsfhigh <- subset_samples(dorms_rare, sex_sheetwashfreq %in% c("male high", "female high"))

##Create new metadata for each condition##
#female low/medium vs female high#
samp_dat_wdiv_f_lowhigh <- subset(samp_dat_wdiv, sex_sheetwashfreq %in% c("female high", "female low"))
#male low/medium vs male high#
samp_dat_wdiv_m_lowhigh <- subset(samp_dat_wdiv, sex_sheetwashfreq %in% c("male high", "male low"))
#female low/medium vs male low/medium#
samp_dat_wdiv_mlowflow <- subset(samp_dat_wdiv, sex_sheetwashfreq %in% c("male low", "female low"))
#female high vs male high#
samp_dat_wdiv_mhighfhigh <- subset(samp_dat_wdiv, sex_sheetwashfreq %in% c("female high", "male high"))

###FEMALE LOW vs FEMALE HIGH###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_f_lowvshigh, method="unifrac")
pcoa_unifrac <- ordinate(phylo_f_lowvshigh, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_f_lowvshigh, pcoa_unifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_f_lowvshigh, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_f_lowhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_F_highvslow_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_f_lowmedvshigh, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_f_lowvshigh, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_f_lowvshigh, pcoa_wunifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_f_lowvshigh, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_f_lowhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_F_highvslow_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_f_lowvshigh, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_f_lowvshigh, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_f_lowvshigh, pcoa_jaccard, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_f_lowvshigh)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwashfreq, data=samp_dat_wdiv_f_lowhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_F_highvslow_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_f_lowvshigh, method="bray")
pcoa_bc <- ordinate(phylo_f_lowvshigh, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_f_lowvshigh, pcoa_bc, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_f_lowvshigh)), method="bray")
adonis2(dm_bray ~ sex_sheetwashfreq, data=samp_dat_wdiv_f_lowhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_F_highvslow_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)



###MALE LOW vs MALE HIGH###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_m_lowvshigh, method="unifrac")
pcoa_unifrac <- ordinate(phylo_m_lowvshigh, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_m_lowvshigh, pcoa_unifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_m_lowvshigh, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_m_lowhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_M_highvslow_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_m_lowvshigh, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_m_lowvshigh, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_m_lowvshigh, pcoa_wunifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_m_lowvshigh, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_m_lowhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_M_highvslow_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_m_lowvshigh, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_m_lowvshigh, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_m_lowvshigh, pcoa_jaccard, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_m_lowvshigh)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwashfreq, data=samp_dat_wdiv_m_lowhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_M_highvslow_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_m_lowvshigh, method="bray")
pcoa_bc <- ordinate(phylo_m_lowvshigh, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_m_lowvshigh, pcoa_bc, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_m_lowvshigh)), method="bray")
adonis2(dm_bray ~ sex_sheetwashfreq, data=samp_dat_wdiv_m_lowhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_M_highvslow_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)



###MALE LOW vs FEMALE LOW###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_mlowvsflow, method="unifrac")
pcoa_unifrac <- ordinate(phylo_mlowvsflow, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_mlowvsflow, pcoa_unifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_mlowvsflow, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_mlowflow)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_low_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_mlowvsflow, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_mlowvsflow, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_mlowvsflow, pcoa_wunifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_mlowvsflow, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_mlowflow)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_low_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_mlowvsflow, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_mlowvsflow, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_mlowvsflow, pcoa_jaccard, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_mlowvsflow)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwashfreq, data=samp_dat_wdiv_mlowflow)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_low_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_mlowvsflow, method="bray")
pcoa_bc <- ordinate(phylo_mlowvsflow, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_mlowvsflow, pcoa_bc, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_mlowvsflow)), method="bray")
adonis2(dm_bray ~ sex_sheetwashfreq, data=samp_dat_wdiv_mlowflow)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_low_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)




###MALE HIGH vs FEMALE HIGH###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_mhighvsfhigh, method="unifrac")
pcoa_unifrac <- ordinate(phylo_mhighvsfhigh, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_mhighvsfhigh, pcoa_unifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_mhighvsfhigh, weighted=FALSE)
adonis2(unifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_mhighfhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_high_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_mhighvsfhigh, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_mhighvsfhigh, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_mhighvsfhigh, pcoa_wunifrac, color = "sex_sheetwashfreq") +
  labs(col="Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_mhighvsfhigh, weighted=TRUE)
adonis2(wunifrac_dm ~ sex_sheetwashfreq, data=samp_dat_wdiv_mhighfhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_high_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_mhighvsfhigh, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_mhighvsfhigh, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_mhighvsfhigh, pcoa_jaccard, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_mhighvsfhigh)), method="jaccard")
adonis2(dm_jaccard ~ sex_sheetwashfreq, data=samp_dat_wdiv_mhighfhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_high_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_mhighvsfhigh, method="bray")
pcoa_bc <- ordinate(phylo_mhighvsfhigh, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_mhighvsfhigh, pcoa_bc, color = "sex_sheetwashfreq") + labs(col = "Sheet Wash Frequency Gender Groups") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_mhighvsfhigh)), method="bray")
adonis2(dm_bray ~ sex_sheetwashfreq, data=samp_dat_wdiv_mhighfhigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_2B/AIM_2B_PT5/2B_5_high_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)
