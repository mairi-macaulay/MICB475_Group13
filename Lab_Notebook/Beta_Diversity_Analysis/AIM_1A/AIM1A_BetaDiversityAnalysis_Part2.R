####AIM 1A Beta Diversity Analysis####

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

##Create new phyloseq for each condition##
#no medium#
phylo_nomedium <- subset_samples(dorms_rare, sheetwashfreq_binned %in% c("low", "high"))
#no low#
phylo_nolow <- subset_samples(dorms_rare, sheetwashfreq_binned %in% c("medium", "high"))
#no high#
phylo_nohigh <- subset_samples(dorms_rare, sheetwashfreq_binned %in% c("low", "medium"))

##Create new metadata for each condition##
#no medium#
samp_dat_wdiv_nomed <- subset(samp_dat_wdiv, sheetwashfreq_binned !="medium")
#no low#
samp_dat_wdiv_nolow <- subset(samp_dat_wdiv, sheetwashfreq_binned !="low")
#no high#
samp_dat_wdiv_nohigh <- subset(samp_dat_wdiv, sheetwashfreq_binned !="high")


###LOW vs HIGH###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_nomedium, method="unifrac")
pcoa_unifrac <- ordinate(phylo_nomedium, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_nomedium, pcoa_unifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_nomedium, weighted=FALSE)
adonis2(unifrac_dm ~ sheetwashfreq_binned, data=samp_dat_wdiv_nomed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2LH_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_nomedium, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_nomedium, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_nomedium, pcoa_wunifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_nomedium, weighted=TRUE)
adonis2(wunifrac_dm ~ sheetwashfreq_binned, data=samp_dat_wdiv_nomed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2LH_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_nomedium, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_nomedium, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_nomedium, pcoa_jaccard, color = "sheetwashfreq_binned") + labs(col = "Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_nomedium)), method="jaccard")
adonis2(dm_jaccard ~ sheetwashfreq_binned, data=samp_dat_wdiv_nomed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2LH_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_nomedium, method="bray")
pcoa_bc <- ordinate(phylo_nomedium, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_nomedium, pcoa_bc, color = "sheetwashfreq_binned") + labs(col = "Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_nomedium)), method="bray")
adonis2(dm_bray ~ sheetwashfreq_binned, data=samp_dat_wdiv_nomed)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2LH_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)


###MEDIUM vs HIGH###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_nolow, method="unifrac")
pcoa_unifrac <- ordinate(phylo_nolow, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_nolow, pcoa_unifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_nolow, weighted=FALSE)
adonis2(unifrac_dm ~ sheetwashfreq_binned, data=samp_dat_wdiv_nolow)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2MH_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_nolow, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_nolow, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_nolow, pcoa_wunifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_nolow, weighted=TRUE)
adonis2(wunifrac_dm ~ sheetwashfreq_binned, data=samp_dat_wdiv_nolow)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2MH_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_nolow, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_nolow, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_nolow, pcoa_jaccard, color = "sheetwashfreq_binned") + labs(col = "Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_nolow)), method="jaccard")
adonis2(dm_jaccard ~ sheetwashfreq_binned, data=samp_dat_wdiv_nolow)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2MH_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_nolow, method="bray")
pcoa_bc <- ordinate(phylo_nolow, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_nolow, pcoa_bc, color = "sheetwashfreq_binned") + labs(col = "Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_nolow)), method="bray")
adonis2(dm_bray ~ sheetwashfreq_binned, data=samp_dat_wdiv_nolow)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2MH_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)



###LOW vs MEDIUM###

##Unweighted Unifrac##
unifrac_dm <- distance(phylo_nohigh, method="unifrac")
pcoa_unifrac <- ordinate(phylo_nohigh, method="PCoA", distance=unifrac_dm)
gg_pcoa_unifrac <- plot_ordination(phylo_nohigh, pcoa_unifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_unifrac

#Unweighted Unifrac PERMANOVA Test#
set.seed(1) 
unifrac_dm <- UniFrac(phylo_nohigh, weighted=FALSE)
adonis2(unifrac_dm ~ sheetwashfreq_binned, data=samp_dat_wdiv_nohigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2LM_unweighted_pcoa.png"
       , gg_pcoa_unifrac
       , height=4, width=5)


##Weighted Unifrac##
wu_dm <- distance(phylo_nohigh, method="wunifrac")
pcoa_wunifrac <- ordinate(phylo_nohigh, method="PCoA", distance=wu_dm)
gg_pcoa_wunifrac <- plot_ordination(phylo_nohigh, pcoa_wunifrac, color = "sheetwashfreq_binned") +
  labs(col="Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_wunifrac

#Weighted Unifrac PERMANOVA Test#
set.seed(1) 
wunifrac_dm <- UniFrac(phylo_nohigh, weighted=TRUE)
adonis2(wunifrac_dm ~ sheetwashfreq_binned, data=samp_dat_wdiv_nohigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2LM_weighted_pcoa.png"
       , gg_pcoa_wunifrac
       , height=4, width=5)


##Jaccard##
j_dm <- distance(phylo_nohigh, method = "jaccard", binary = TRUE)
pcoa_jaccard <- ordinate(phylo_nohigh, method="PCoA", distance=j_dm)
gg_pcoa_jaccard <- plot_ordination(phylo_nohigh, pcoa_jaccard, color = "sheetwashfreq_binned") + labs(col = "Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_jaccard

#Jaccard Diversity PERMANOVA Test#
set.seed(1) 
dm_jaccard <- vegdist(t(otu_table(phylo_nohigh)), method="jaccard")
adonis2(dm_jaccard ~ sheetwashfreq_binned, data=samp_dat_wdiv_nohigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2LM_jaccard_pcoa.png"
       , gg_pcoa_jaccard
       , height=4, width=5)


##Bray Curtis##
bc_dm <- distance(phylo_nohigh, method="bray")
pcoa_bc <- ordinate(phylo_nohigh, method="PCoA", distance=bc_dm)
gg_pcoa_bc <- plot_ordination(phylo_nohigh, pcoa_bc, color = "sheetwashfreq_binned") + labs(col = "Sheet Wash Frequency") + stat_ellipse()
gg_pcoa_bc

#Bray Curtis Diversity PERMANOVA Test#
set.seed(1) 
dm_bray <- vegdist(t(otu_table(phylo_nohigh)), method="bray")
adonis2(dm_bray ~ sheetwashfreq_binned, data=samp_dat_wdiv_nohigh)

#save PCoA#
ggsave(filename = "Lab_Notebook/Beta_Diversity_Analysis/AIM_1A/1A_2LM_bc_pcoa.png"
       , gg_pcoa_bc
       , height=4, width=5)

