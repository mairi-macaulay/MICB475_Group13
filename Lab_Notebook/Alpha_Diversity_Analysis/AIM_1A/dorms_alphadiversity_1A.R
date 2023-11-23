####ALPHA DIVERSITY (1A - Sheet Washing Frequency)####

#set seed
set.seed(1)

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
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_shannon_sheetwashfreq <- ggplot(samp_dat_wdiv, aes(x=`sheetwashfreq_binned`, y=Shannon)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency")
gg_shannon_sheetwashfreq

#One-way ANOVA
# Set up our linear model 
lm_shannon_vs_sheetwashfreq <- lm(Shannon ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv)
# Calculate AOV
anova_shannon_vs_sheetwashfreq <- aov(lm_shannon_vs_sheetwashfreq)
# Summarize to determine if there are significant differences
summary(anova_shannon_vs_sheetwashfreq)
# Determine which groups are significant
tukey_sum_shannon <- TukeyHSD(anova_shannon_vs_sheetwashfreq)

#Results: no significance found 

# save plot file 
ggsave(filename = "sheetwashfreq_plot_shannon.png"
       , gg_shannon_sheetwashfreq
       , height=4, width=6)

#####Observed##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_observed_sheetwashfreq <- ggplot(samp_dat_wdiv, aes(x=`sheetwashfreq_binned`, y=Observed)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency")
gg_observed_sheetwashfreq

#ANOVA
# Set up our linear model 
lm_observed_vs_sheetwashfreq <- lm(Observed ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv)
# Calculate AOV
anova_observed_vs_sheetwashfreq <- aov(lm_observed_vs_sheetwashfreq)
# Summarize to determine if there are significant differences
summary(anova_observed_vs_sheetwashfreq)
# Determine which groups are significant
tukey_sum_observed <- TukeyHSD(anova_observed_vs_sheetwashfreq)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_plot_observed.png"
       , gg_observed_sheetwashfreq
       , height=4, width=6)

#####Chao1##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_chao1_sheetwashfreq <- ggplot(samp_dat_wdiv, aes(x=`sheetwashfreq_binned`, y=Chao1)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency")
gg_chao1_sheetwashfreq

#ANOVA
# Set up our linear model 
lm_chao1_vs_sheetwashfreq <- lm(Chao1 ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv)
# Calculate AOV
anova_chao1_vs_sheetwashfreq <- aov(lm_chao1_vs_sheetwashfreq)
# Summarize to determine if there are significant differences
summary(anova_chao1_vs_sheetwashfreq)
# Determine which groups are significant
tukey_sum_chao1 <- TukeyHSD(anova_chao1_vs_sheetwashfreq)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_plot_chao1.png"
       , gg_chao1_sheetwashfreq
       , height=4, width=6)

#####ACE##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_ace_sheetwashfreq <- ggplot(samp_dat_wdiv, aes(x=`sheetwashfreq_binned`, y=ACE)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency")
gg_ace_sheetwashfreq

#ANOVA
# Set up our linear model 
lm_ace_vs_sheetwashfreq <- lm(ACE ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv)
# Calculate AOV
anova_ace_vs_sheetwashfreq <- aov(lm_ace_vs_sheetwashfreq)
# Summarize to determine if there are significant differences
summary(anova_ace_vs_sheetwashfreq)
# Determine which groups are significant
tukey_sum_ace <- TukeyHSD(anova_ace_vs_sheetwashfreq)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_plot_ace.png"
       , gg_ace_sheetwashfreq
       , height=4, width=6)

#####Simpson##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_simpson_sheetwashfreq <- ggplot(samp_dat_wdiv, aes(x=`sheetwashfreq_binned`, y=Simpson)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency")
gg_simpson_sheetwashfreq

#ANOVA
# Set up our linear model 
lm_simpson_vs_sheetwashfreq <- lm(Simpson ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv)
# Calculate AOV
anova_simpson_vs_sheetwashfreq <- aov(lm_simpson_vs_sheetwashfreq)
# Summarize to determine if there are significant differences
summary(anova_simpson_vs_sheetwashfreq)
# Determine which groups are significant
tukey_sum_simpson <- TukeyHSD(anova_simpson_vs_sheetwashfreq)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_plot_simpson.png"
       , gg_simpson_sheetwashfreq
       , height=4, width=6)

#####InvSimpson##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_invsimpson_sheetwashfreq <- ggplot(samp_dat_wdiv, aes(x=`sheetwashfreq_binned`, y=InvSimpson)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency")
gg_invsimpson_sheetwashfreq

#ANOVA
# Set up our linear model 
lm_invsimp_vs_sheetwashfreq <- lm(InvSimpson ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv)
# Calculate AOV
anova_invsimp_vs_sheetwashfreq <- aov(lm_invsimp_vs_sheetwashfreq)
# Summarize to determine if there are significant differences
summary(anova_invsimp_vs_sheetwashfreq)
# Determine which groups are significant
tukey_sum_invsimp <- TukeyHSD(anova_invsimp_vs_sheetwashfreq)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_plot_invsimpson.png"
       , gg_invsimpson_sheetwashfreq
       , height=4, width=6)

#####Fisher##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_fisher_sheetwashfreq <- ggplot(samp_dat_wdiv, aes(x=`sheetwashfreq_binned`, y=Fisher)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency")
gg_fisher_sheetwashfreq

#ANOVA
# Set up our linear model 
lm_fisher_vs_sheetwashfreq <- lm(Fisher ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv)
# Calculate AOV
anova_fisher_vs_sheetwashfreq <- aov(lm_fisher_vs_sheetwashfreq)
# Summarize to determine if there are significant differences
summary(anova_fisher_vs_sheetwashfreq)
# Determine which groups are significant
tukey_sum_fisher <- TukeyHSD(anova_fisher_vs_sheetwashfreq)

#Results: no significance found

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
plot.pd_sheetwashfreq <- ggplot(sample_data(dorms_rare), aes(sheetwashfreq_binned, PD)) +
  geom_boxplot() +
  xlab("Sheet Washing Frequency") +
  ylab("Phylogenetic Diversity")
# view plot
plot.pd_sheetwashfreq

# save plot file 
ggsave(filename = "sheetwashfreq_plot_pd.png"
       , plot.pd_sheetwashfreq
       , height=4, width=6)

