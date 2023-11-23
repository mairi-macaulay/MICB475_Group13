####ALPHA DIVERSITY (1B - Shower Recency)####

#set seed
set.seed(1)

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)

##Load in RData
load("../../Phyloseq/dorms_rare_showerrecency.RData")

#### Alpha diversity ######
#view all metrics
plot_richness(dorms_rare)

#extract information and combine with alpha diversity metrics
alphadiv <- estimate_richness(dorms_rare)
samp_dat <- sample_data(dorms_rare) #extract out metadata
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

#####Shannon#####
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_shannon_showerrecency <- ggplot(samp_dat_wdiv, aes(x=`last_shower_binned`, y=Shannon)) +
  geom_boxplot() +
  xlab("Shower Recency")
gg_shannon_showerrecency

#T-test
t.test(samp_dat_wdiv$Shannon ~ samp_dat_wdiv$last_shower_binned)

#Results: p-value = 0.8094 > 0.05 = not significant

# save plot file 
ggsave(filename = "showerrecency_plot_shannon.png"
       , gg_shannon_showerrecency
       , height=4, width=6)

#####Observed##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_observed_showerrecency <- ggplot(samp_dat_wdiv, aes(x=`last_shower_binned`, y=Observed)) +
  geom_boxplot() +
  xlab("Shower Recency")
gg_observed_showerrecency

#T-test
t.test(samp_dat_wdiv$Observed ~ samp_dat_wdiv$last_shower_binned)

#Results: p-value = 0.4466 > 0.05 = not significant

# save plot file 
ggsave(filename = "showerrecency_plot_observed.png"
       , gg_observed_showerrecency
       , height=4, width=6)

#####Chao1##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_chao1_showerrecency <- ggplot(samp_dat_wdiv, aes(x=`last_shower_binned`, y=Chao1)) +
  geom_boxplot() +
  xlab("Shower Recency")
gg_chao1_showerrecency

#T-test
t.test(samp_dat_wdiv$Chao1 ~ samp_dat_wdiv$last_shower_binned)

#Results: p-value = 0.4827 > 0.05 = not significant

# save plot file 
ggsave(filename = "showerrecency_plot_chao1.png"
       , gg_chao1_showerrecency
       , height=4, width=6)

#####ACE##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_ace_showerrecency <- ggplot(samp_dat_wdiv, aes(x=`last_shower_binned`, y=ACE)) +
  geom_boxplot() +
  xlab("Shower Recency")
gg_ace_showerrecency

#T-test
t.test(samp_dat_wdiv$ACE ~ samp_dat_wdiv$last_shower_binned)

#Results: p-value = 0.4849 > 0.05 = not significant

# save plot file 
ggsave(filename = "showerrecency_plot_ace.png"
       , gg_ace_showerrecency
       , height=4, width=6)

#####Simpson##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_simpson_showerrecency <- ggplot(samp_dat_wdiv, aes(x=`last_shower_binned`, y=Simpson)) +
  geom_boxplot() +
  xlab("Shower Recency")
gg_simpson_showerrecency

#T-test
t.test(samp_dat_wdiv$Simpson ~ samp_dat_wdiv$last_shower_binned)

#Results: p-value = 0.5755 > 0.05 = not significant

# save plot file 
ggsave(filename = "showerrecency_plot_simpson.png"
       , gg_simpson_showerrecency
       , height=4, width=6)

#####InvSimpson##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_invsimpson_showerrecency <- ggplot(samp_dat_wdiv, aes(x=`last_shower_binned`, y=InvSimpson)) +
  geom_boxplot() +
  xlab("Shower Recency")
gg_invsimpson_showerrecency

#T-test
t.test(samp_dat_wdiv$InvSimpson ~ samp_dat_wdiv$last_shower_binned)

#Results: p-value = 0.5335 > 0.05 = not significant

# save plot file 
ggsave(filename = "showerrecency_plot_invsimpson.png"
       , gg_invsimpson_showerrecency
       , height=4, width=6)

#####Fisher##### 
# save into object, define x-axis as last_shower_binned and relabel, add a boxplot
gg_fisher_showerrecency <- ggplot(samp_dat_wdiv, aes(x=`last_shower_binned`, y=Fisher)) +
  geom_boxplot() +
  xlab("Shower Recency")
gg_fisher_showerrecency

#T-test
t.test(samp_dat_wdiv$Fisher ~ samp_dat_wdiv$last_shower_binned)

#Results: p-value = 0.4621 > 0.05 = not significant

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

# save plot file 
ggsave(filename = "showerrecency_plot_pd.png"
       , plot.pd_showerrecency
       , height=4, width=6)

