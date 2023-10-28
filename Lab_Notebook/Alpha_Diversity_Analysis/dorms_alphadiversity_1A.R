####ALPHA DIVERSITY (1A - Sheet Washing Frequency)####

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

##Load in RData
load("../Phyloseq/dorms_rare_sheetwashfreq.RData")
load("../Phyloseq/dorms_final_sheetwashfreq.RData")

#### Alpha diversity ######
#view all metrics
plot_richness(dorms_rare)

#extract information
alphadiv <- estimate_richness(dorms_rare)
samp_dat <- sample_data(dorms_rare) #extract out metadata
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

#####Shannon#####
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_shannon_sheetwashfreq <- plot_richness(dorms_rare, x = "sheetwashfreq_binned", measures = c("Shannon")) +
  xlab("Sheet Wash Frequency") +
  geom_boxplot()
gg_shannon_sheetwashfreq

#ANOVA

# save plot file 
ggsave(filename = "sheetwashfreq_plot_shannon.png"
       , gg_shannon_sheetwashfreq
       , height=4, width=6)

#####Observed##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_observed_sheetwashfreq <- plot_richness(dorms_rare, x = "sheetwashfreq_binned", measures = c("Observed")) +
  xlab("Sheet Wash Frequency") +
  geom_boxplot()
gg_observed_sheetwashfreq

#ANOVA

# save plot file 
ggsave(filename = "sheetwashfreq_plot_observed.png"
       , gg_observed_sheetwashfreq
       , height=4, width=6)

#####Chao1##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_chao1_sheetwashfreq <- plot_richness(dorms_rare, x = "sheetwashfreq_binned", measures = c("Chao1")) +
  xlab("Sheet Wash Frequency") +
  geom_boxplot()
gg_chao1_sheetwashfreq

#ANOVA

# save plot file 
ggsave(filename = "sheetwashfreq_plot_chao1.png"
       , gg_chao1_sheetwashfreq
       , height=4, width=6)

#####ACE##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_ace_sheetwashfreq <- plot_richness(dorms_rare, x = "sheetwashfreq_binned", measures = c("ACE")) +
  xlab("Sheet Wash Frequency") +
  geom_boxplot()
gg_ace_sheetwashfreq

#ANOVA

# save plot file 
ggsave(filename = "sheetwashfreq_plot_ace.png"
       , gg_ace_sheetwashfreq
       , height=4, width=6)

#####Simpson##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_simpson_sheetwashfreq <- plot_richness(dorms_rare, x = "sheetwashfreq_binned", measures = c("Simpson")) +
  xlab("Sheet Wash Frequency") +
  geom_boxplot()
gg_simpson_sheetwashfreq

#ANOVA

# save plot file 
ggsave(filename = "sheetwashfreq_plot_simpson.png"
       , gg_simpson_sheetwashfreq
       , height=4, width=6)

#####InvSimpson##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_invsimpson_sheetwashfreq <- plot_richness(dorms_rare, x = "sheetwashfreq_binned", measures = c("InvSimpson")) +
  xlab("Sheet Wash Frequency") +
  geom_boxplot()
gg_invsimpson_sheetwashfreq

#ANOVA

# save plot file 
ggsave(filename = "sheetwashfreq_plot_invsimpson.png"
       , gg_invsimpson_sheetwashfreq
       , height=4, width=6)

#####Fisher##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_fisher_sheetwashfreq <- plot_richness(dorms_rare, x = "sheetwashfreq_binned", measures = c("Fisher")) +
  xlab("Sheet Wash Frequency") +
  geom_boxplot()
gg_fisher_sheetwashfreq

#ANOVA

# save plot file 
ggsave(filename = "sheetwashfreq_plot_fisher.png"
       , gg_fisher_sheetwashfreq
       , height=4, width=6)

#####Phylogenetic diversity#####
# calculate Faith's phylogenetic diversity as PD
phylo_dist_sheetwashfreq <- pd(t(otu_table(dorms_rare)), phy_tree(dorms_rare),
                 include.root=F)
# add PD to metadata table
sample_data(dorms_rare)$PD <- phylo_dist_sheetwashfreq$PD
# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(dorms_rare), aes(sheetwashfreq_binned, PD)) +
  geom_boxplot() +
  xlab("Sheet Washing Frequency") +
  ylab("Phylogenetic Diversity")
# view plot
plot.pd_sheetwashfreq

#ANOVA


# save plot file 
ggsave(filename = "sheetwashfreq_plot_pd.png"
       , plot.pd_sheetwashfreq
       , height=4, width=6)

