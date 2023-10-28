####ALPHA DIVERSITY (1B - Shower Recency)####

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

##Load in RData
load("../Phyloseq/dorms_rare_showerrecency.RData")
load("../Phyloseq/dorms_final_showerrecency.RData")

#### Alpha diversity ######
#view all metrics
plot_richness(dorms_rare)

#extract information
alphadiv <- estimate_richness(dorms_rare)
samp_dat <- sample_data(dorms_rare) #extract out metadata
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

#####Shannon#####
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_shannon_showerrecency <- plot_richness(dorms_rare, x = "last_shower_binned", measures = c("Shannon")) +
  xlab("Shower Recency") +
  geom_boxplot()
gg_shannon_showerrecency

#ANOVA

# save plot file 
ggsave(filename = "showerrecency_plot_shannon.png"
       , gg_shannon_showerrecency
       , height=4, width=6)

#####Observed##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_observed_showerrecency <- plot_richness(dorms_rare, x = "last_shower_binned", measures = c("Observed")) +
  xlab("Shower Recency") +
  geom_boxplot()
gg_observed_showerrecency

#ANOVA

# save plot file 
ggsave(filename = "showerrecency_plot_observed.png"
       , gg_observed_showerrecency
       , height=4, width=6)

#####Chao1##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_chao1_showerrecency <- plot_richness(dorms_rare, x = "last_shower_binned", measures = c("Chao1")) +
  xlab("Shower Recency") +
  geom_boxplot()
gg_chao1_showerrecency

#ANOVA

# save plot file 
ggsave(filename = "showerrecency_plot_chao1.png"
       , gg_chao1_showerrecency
       , height=4, width=6)

#####ACE##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_ace_showerrecency <- plot_richness(dorms_rare, x = "last_shower_binned", measures = c("ACE")) +
  xlab("Shower Recency") +
  geom_boxplot()
gg_ace_showerrecency

#ANOVA

# save plot file 
ggsave(filename = "showerrecency_plot_ace.png"
       , gg_ace_showerrecency
       , height=4, width=6)

#####Simpson##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_simpson_showerrecency <- plot_richness(dorms_rare, x = "last_shower_binned", measures = c("Simpson")) +
  xlab("Shower Recency") +
  geom_boxplot()
gg_simpson_showerrecency

#ANOVA

# save plot file 
ggsave(filename = "showerrecency_plot_simpson.png"
       , gg_simpson_showerrecency
       , height=4, width=6)

#####InvSimpson##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_invsimpson_showerrecency <- plot_richness(dorms_rare, x = "last_shower_binned", measures = c("InvSimpson")) +
  xlab("Shower Recency") +
  geom_boxplot()
gg_invsimpson_showerrecency

#ANOVA

# save plot file 
ggsave(filename = "showerrecency_plot_invsimpson.png"
       , gg_invsimpson_showerrecency
       , height=4, width=6)

#####Fisher##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_fisher_showerrecency <- plot_richness(dorms_rare, x = "last_shower_binned", measures = c("Fisher")) +
  xlab("Shower Recency") +
  geom_boxplot()
gg_fisher_showerrecency

#ANOVA

# save plot file 
ggsave(filename = "showerrecency_plot_fisher.png"
       , gg_fisher_showerrecency
       , height=4, width=6)

#####Phylogenetic diversity#####
# calculate Faith's phylogenetic diversity as PD
phylo_dist_showerrecency <- pd(t(otu_table(dorms_rare)), phy_tree(dorms_rare),
                 include.root=F)
# add PD to metadata table
sample_data(dorms_rare)$PD <- phylo_dist_showerrecency$PD
# plot any metadata category against the PD
plot.pd_showerrecency <- ggplot(sample_data(dorms_rare), aes(last_shower_binned, PD)) +
  geom_boxplot() +
  xlab("Shower Recency") +
  ylab("Phylogenetic Diversity")
# view plot
plot.pd_showerrecency

#ANOVA


# save plot file 
ggsave(filename = "showerrecency_plot_pd.png"
       , plot.pd_showerrecency
       , height=4, width=6)

