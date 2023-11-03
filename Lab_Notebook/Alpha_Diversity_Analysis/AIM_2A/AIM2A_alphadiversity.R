####ALPHA DIVERSITY (2A - GENDER)####

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)

##Load in RData
load("../../Phyloseq/dorms_rare_sheetwashfreq.RData")

#### Alpha diversity ######
#view all metrics
plot_richness(dorms_rare)

#extract information and combine with alpha diversity metrics
alphadiv <- estimate_richness(dorms_rare)
samp_dat <- sample_data(dorms_rare) #extract out metadata
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

#####Shannon#####
# save into object, define x-axis as sex and relabel, add a boxplot
gg_shannon_gender <- ggplot(samp_dat_wdiv, aes(x=`sex`, y=Shannon)) +
  geom_boxplot() +
  xlab("Gender")
gg_shannon_gender

#T-test
t.test(samp_dat_wdiv$Shannon ~ samp_dat_wdiv$sex)
#p-value = 0.5554

# save plot file 
ggsave(filename = "2A_plot_shannon.png"
       , gg_shannon_gender
       , height=4, width=6)

#####Observed##### 
# save into object, define x-axis as sex and relabel, add a boxplot
gg_observed_gender <- ggplot(samp_dat_wdiv, aes(x=`sex`, y=Observed)) +
  geom_boxplot() +
  xlab("Gender")
gg_observed_gender

#T-test
t.test(samp_dat_wdiv$Observed ~ samp_dat_wdiv$sex)
#p-value = 0.7671

# save plot file 
ggsave(filename = "2A_plot_observed.png"
       , gg_observed_gender
       , height=4, width=6)

#####Chao1##### 
# save into object, define x-axis as sex and relabel, add a boxplot
gg_chao1_gender <- ggplot(samp_dat_wdiv, aes(x=`sex`, y=Chao1)) +
  geom_boxplot() +
  xlab("Gender")
gg_chao1_gender

#T-test
t.test(samp_dat_wdiv$Chao1 ~ samp_dat_wdiv$sex)
#p-value = 0.7531

# save plot file 
ggsave(filename = "2A_plot_chao1.png"
       , gg_chao1_gender
       , height=4, width=6)

#####ACE##### 
# save into object, define x-axis as sex and relabel, add a boxplot
gg_ace_gender <- ggplot(samp_dat_wdiv, aes(x=`sex`, y=ACE)) +
  geom_boxplot() +
  xlab("Gender")
gg_ace_gender

#T-test
t.test(samp_dat_wdiv$ACE ~ samp_dat_wdiv$sex)
#p-value = 0.7161

# save plot file 
ggsave(filename = "2A_plot_ace.png"
       , gg_ace_gender
       , height=4, width=6)

#####Simpson##### 
# save into object, define x-axis as sex and relabel, add a boxplot
gg_simpson_gender <- ggplot(samp_dat_wdiv, aes(x=`sex`, y=Simpson)) +
  geom_boxplot() +
  xlab("Gender")
gg_simpson_gender

#T-test
t.test(samp_dat_wdiv$Simpson ~ samp_dat_wdiv$sex)
#p-value = 0.1682

# save plot file 
ggsave(filename = "2A_plot_simpson.png"
       , gg_simpson_gender
       , height=4, width=6)

#####InvSimpson##### 
# save into object, define x-axis as sex and relabel, add a boxplot
gg_invsimpson_gender <- ggplot(samp_dat_wdiv, aes(x=`sex`, y=InvSimpson)) +
  geom_boxplot() +
  xlab("Gender")
gg_invsimpson_gender

#T-test
t.test(samp_dat_wdiv$InvSimpson ~ samp_dat_wdiv$sex)
#p-value = 0.636

# save plot file 
ggsave(filename = "2A_plot_invsimpson.png"
       , gg_invsimpson_gender
       , height=4, width=6)

#####Fisher##### 
# save into object, define x-axis as sex and relabel, add a boxplot
gg_fisher_gender <- ggplot(samp_dat_wdiv, aes(x=`sex`, y=Fisher)) +
  geom_boxplot() +
  xlab("Gender")
gg_fisher_gender

#T-test
t.test(samp_dat_wdiv$Fisher ~ samp_dat_wdiv$sex)
#p-value = 0.6716

# save plot file 
ggsave(filename = "2A_plot_fisher.png"
       , gg_fisher_gender
       , height=4, width=6)

#####Phylogenetic diversity#####
# calculate Faith's phylogenetic diversity as PD
phylo_dist_gender <- pd(t(otu_table(dorms_rare)), phy_tree(dorms_rare),
                               include.root=F)
# add PD to metadata table
sample_data(dorms_rare)$PD <- phylo_dist_gender$PD
# plot any metadata category against the PD
plot.pd_gender <- ggplot(sample_data(dorms_rare), aes(sex, PD)) +
  geom_boxplot() +
  xlab("Gender") +
  ylab("Phylogenetic Diversity")
# view plot
plot.pd_gender

# save plot file 
ggsave(filename = "2A_plot_pd.png"
       , plot.pd_gender
       , height=4, width=6)
