####ALPHA DIVERSITY (2B - Sheet Washing Frequency and Gender)####

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)

##Load in RData
load("../../Phyloseq/AIM_2B_phyloseq/dorms_rare_sheetwashfreq_female.RData")
load("../../Phyloseq/AIM_2B_phyloseq/dorms_rare_sheetwashfreq_male.RData")

#### Alpha diversity - FEMALE######
#view all metrics for female
plot_richness(dorms_rare_sheetwashfreq_female)

#extract information and combine with alpha diversity metrics
alphadiv_female <- estimate_richness(dorms_rare_sheetwashfreq_female)
samp_dat_female <- sample_data(dorms_rare_sheetwashfreq_female) #extract out metadata
samp_dat_wdiv_female <- data.frame(samp_dat_female, alphadiv_female)

#####Shannon#####
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_shannon_sheetwashfreq_female <- ggplot(samp_dat_wdiv_female, aes(x=`sheetwashfreq_binned`, y=Shannon)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Shannon Metrics for Female Sheet Wash Frequency")
gg_shannon_sheetwashfreq_female

#One-way ANOVA
# Set up our linear model 
lm_shannon_vs_sheetwashfreq_female <- lm(Shannon ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_female)
# Calculate AOV
anova_shannon_vs_sheetwashfreq_female <- aov(lm_shannon_vs_sheetwashfreq_female)
# Summarize to determine if there are significant differences
summary(anova_shannon_vs_sheetwashfreq_female)
# Determine which groups are significant
tukey_sum_shannon_female <- TukeyHSD(anova_shannon_vs_sheetwashfreq_female)

# Kruskal-wallis test
kruskal_shannon_female <- kruskal.test(Shannon ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_female)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_shannon_vs_sheetwashfreq_female_log <- lm(log(Shannon) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_female)
anova_shannon_vs_sheetwashfreq_female_log <- aov(lm_shannon_vs_sheetwashfreq_female_log)
summary(anova_shannon_vs_sheetwashfreq_female_log)
TukeyHSD(anova_shannon_vs_sheetwashfreq_female_log)

#Results: no significance

# save plot file 
ggsave(filename = "sheetwashfreq_female_plot_shannon.png"
       , gg_shannon_sheetwashfreq_female
       , height=4, width=6)

#####Observed##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_observed_sheetwashfreq_female <- ggplot(samp_dat_wdiv_female, aes(x=`sheetwashfreq_binned`, y=Observed)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Observed Metrics for Female Sheet Wash Frequency")
gg_observed_sheetwashfreq_female

#ANOVA
# Set up our linear model 
lm_observed_vs_sheetwashfreq_female <- lm(Observed ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_female)
# Calculate AOV
anova_observed_vs_sheetwashfreq_female <- aov(lm_observed_vs_sheetwashfreq_female)
# Summarize to determine if there are significant differences
summary(anova_observed_vs_sheetwashfreq_female)
# Determine which groups are significant
tukey_sum_observed_female <- TukeyHSD(anova_observed_vs_sheetwashfreq_female)

# Kruskal-wallis test
kruskal_observed_female <- kruskal.test(Observed ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_female)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_observed_vs_sheetwashfreq_female_log <- lm(log(Observed) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_female)
anova_observed_vs_sheetwashfreq_female_log <- aov(lm_observed_vs_sheetwashfreq_female_log)
summary(anova_observed_vs_sheetwashfreq_female_log)
TukeyHSD(anova_observed_vs_sheetwashfreq_female_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_female_plot_observed.png"
       , gg_observed_sheetwashfreq_female
       , height=4, width=6)

#####Chao1##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_chao1_sheetwashfreq_female <- ggplot(samp_dat_wdiv_female, aes(x=`sheetwashfreq_binned`, y=Chao1)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Chao1 Metrics for Female Sheet Wash Frequency")
gg_chao1_sheetwashfreq_female

#ANOVA
# Set up our linear model 
lm_chao1_vs_sheetwashfreq_female <- lm(Chao1 ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_female)
# Calculate AOV
anova_chao1_vs_sheetwashfreq_female <- aov(lm_chao1_vs_sheetwashfreq_female)
# Summarize to determine if there are significant differences
summary(anova_chao1_vs_sheetwashfreq_female)
# Determine which groups are significant
tukey_sum_chao1_female <- TukeyHSD(anova_chao1_vs_sheetwashfreq_female)

# Kruskal-wallis test
kruskal_chao1_female <- kruskal.test(Chao1 ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_female)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_chao1_vs_sheetwashfreq_female_log <- lm(log(Chao1) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_female)
anova_chao1_vs_sheetwashfreq_female_log <- aov(lm_chao1_vs_sheetwashfreq_female_log)
summary(anova_chao1_vs_sheetwashfreq_female_log)
TukeyHSD(anova_chao1_vs_sheetwashfreq_female_log)

# save plot file 
ggsave(filename = "sheetwashfreq_female_plot_chao1.png"
       , gg_chao1_sheetwashfreq_female
       , height=4, width=6)

#####ACE##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_ace_sheetwashfreq_female <- ggplot(samp_dat_wdiv_female, aes(x=`sheetwashfreq_binned`, y=ACE)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "ACE Metrics for Female Sheet Wash Frequency")
gg_ace_sheetwashfreq_female

#ANOVA
# Set up our linear model 
lm_ace_vs_sheetwashfreq_female <- lm(ACE ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_female)
# Calculate AOV
anova_ace_vs_sheetwashfreq_female <- aov(lm_ace_vs_sheetwashfreq_female)
# Summarize to determine if there are significant differences
summary(anova_ace_vs_sheetwashfreq_female)
# Determine which groups are significant
tukey_sum_ace_female <- TukeyHSD(anova_ace_vs_sheetwashfreq_female)

# Kruskal-wallis test
kruskal_ace_female <- kruskal.test(ACE ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_female)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_ace_vs_sheetwashfreq_female_log <- lm(log(ACE) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_female)
anova_ace_vs_sheetwashfreq_female_log <- aov(lm_ace_vs_sheetwashfreq_female_log)
summary(anova_ace_vs_sheetwashfreq_female_log)
TukeyHSD(anova_ace_vs_sheetwashfreq_female_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_female_plot_ace.png"
       , gg_ace_sheetwashfreq_female
       , height=4, width=6)

#####Simpson##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_simpson_sheetwashfreq_female <- ggplot(samp_dat_wdiv_female, aes(x=`sheetwashfreq_binned`, y=Simpson)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Simpson Metrics for Female Sheet Wash Frequency")
gg_simpson_sheetwashfreq_female

#ANOVA
# Set up our linear model 
lm_simpson_vs_sheetwashfreq_female <- lm(Simpson ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_female)
# Calculate AOV
anova_simpson_vs_sheetwashfreq_female <- aov(lm_simpson_vs_sheetwashfreq_female)
# Summarize to determine if there are significant differences
summary(anova_simpson_vs_sheetwashfreq_female)
# Determine which groups are significant
tukey_sum_simpson_female <- TukeyHSD(anova_simpson_vs_sheetwashfreq_female)

# Kruskal-wallis test
kruskal_simpson_female <- kruskal.test(Simpson ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_female)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_simpson_vs_sheetwashfreq_female_log <- lm(log(Simpson) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_female)
anova_simpson_vs_sheetwashfreq_female_log <- aov(lm_simpson_vs_sheetwashfreq_female_log)
summary(anova_simpson_vs_sheetwashfreq_female_log)
TukeyHSD(anova_simpson_vs_sheetwashfreq_female_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_female_plot_simpson.png"
       , gg_simpson_sheetwashfreq_female
       , height=4, width=6)

#####InvSimpson##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_invsimpson_sheetwashfreq_female <- ggplot(samp_dat_wdiv_female, aes(x=`sheetwashfreq_binned`, y=InvSimpson)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "InvSimpson Metrics for Female Sheet Wash Frequency")
gg_invsimpson_sheetwashfreq_female

#ANOVA
# Set up our linear model 
lm_invsimp_vs_sheetwashfreq_female <- lm(InvSimpson ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_female)
# Calculate AOV
anova_invsimp_vs_sheetwashfreq_female <- aov(lm_invsimp_vs_sheetwashfreq_female)
# Summarize to determine if there are significant differences
summary(anova_invsimp_vs_sheetwashfreq_female)
# Determine which groups are significant
tukey_sum_invsimp_female <- TukeyHSD(anova_invsimp_vs_sheetwashfreq_female)

# Kruskal-wallis test
kruskal_invsimp_female <- kruskal.test(InvSimpson ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_female)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_invsimp_vs_sheetwashfreq_female_log <- lm(log(InvSimpson) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_female)
anova_invsimp_vs_sheetwashfreq_female_log <- aov(lm_invsimp_vs_sheetwashfreq_female_log)
summary(anova_invsimp_vs_sheetwashfreq_female_log)
TukeyHSD(anova_invsimp_vs_sheetwashfreq_female_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_female_plot_invsimpson.png"
       , gg_invsimpson_sheetwashfreq_female
       , height=4, width=6)

#####Fisher##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_fisher_sheetwashfreq_female <- ggplot(samp_dat_wdiv_female, aes(x=`sheetwashfreq_binned`, y=Fisher)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Fisher Metrics for Female Sheet Wash Frequency")
gg_fisher_sheetwashfreq_female

#ANOVA
# Set up our linear model 
lm_fisher_vs_sheetwashfreq_female <- lm(Fisher ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_female)
# Calculate AOV
anova_fisher_vs_sheetwashfreq_female <- aov(lm_fisher_vs_sheetwashfreq_female)
# Summarize to determine if there are significant differences
summary(anova_fisher_vs_sheetwashfreq_female)
# Determine which groups are significant
tukey_sum_fisher_female <- TukeyHSD(anova_fisher_vs_sheetwashfreq_female)

# Kruskal-wallis test
kruskal_fisher_female <- kruskal.test(Fisher ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_female)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_fisher_vs_sheetwashfreq_female_log <- lm(log(Fisher) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_female)
anova_fisher_vs_sheetwashfreq_female_log <- aov(lm_fisher_vs_sheetwashfreq_female_log)
summary(anova_fisher_vs_sheetwashfreq_female_log)
TukeyHSD(anova_fisher_vs_sheetwashfreq_female_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_female_plot_fisher.png"
       , gg_fisher_sheetwashfreq_female
       , height=4, width=6)

#####Phylogenetic diversity#####
# calculate Faith's phylogenetic diversity as PD
phylo_dist_sheetwashfreq_female <- pd(t(otu_table(dorms_rare_sheetwashfreq_female)), phy_tree(dorms_rare_sheetwashfreq_female),
                               include.root=F)
# add PD to metadata table
sample_data(dorms_rare_sheetwashfreq_female)$PD <- phylo_dist_sheetwashfreq_female$PD
# plot any metadata category against the PD
plot.pd_sheetwashfreq_female <- ggplot(sample_data(dorms_rare_sheetwashfreq_female), aes(sheetwashfreq_binned, PD)) +
  geom_boxplot() +
  xlab("Sheet Washing Frequency") +
  ylab("Phylogenetic Diversity") +
  labs(title = "Phylogenetic Diversity for Female Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "medium", "high"))
# view plot
plot.pd_sheetwashfreq_female

# save plot file 
ggsave(filename = "sheetwashfreq_female_plot_pd.png"
       , plot.pd_sheetwashfreq_female
       , height=4, width=6)

#### Alpha diversity - MALE######
#view all metrics for male
plot_richness(dorms_rare_sheetwashfreq_male)

#extract information and combine with alpha diversity metrics
alphadiv_male <- estimate_richness(dorms_rare_sheetwashfreq_male)
samp_dat_male <- sample_data(dorms_rare_sheetwashfreq_male) #extract out metadata
samp_dat_wdiv_male <- data.frame(samp_dat_male, alphadiv_male)

#####Shannon#####
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_shannon_sheetwashfreq_male <- ggplot(samp_dat_wdiv_male, aes(x=`sheetwashfreq_binned`, y=Shannon)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Shannon Metrics for Male Sheet Wash Frequency")
gg_shannon_sheetwashfreq_male

#One-way ANOVA
# Set up our linear model 
lm_shannon_vs_sheetwashfreq_male <- lm(Shannon ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_male)
# Calculate AOV
anova_shannon_vs_sheetwashfreq_male <- aov(lm_shannon_vs_sheetwashfreq_male)
# Summarize to determine if there are significant differences
summary(anova_shannon_vs_sheetwashfreq_male)
# Determine which groups are significant
tukey_sum_shannon_male <- TukeyHSD(anova_shannon_vs_sheetwashfreq_male)

# Kruskal-wallis test
kruskal_shannon_male <- kruskal.test(Shannon ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_male)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_shannon_vs_sheetwashfreq_male_log <- lm(log(Shannon) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_male)
anova_shannon_vs_sheetwashfreq_male_log <- aov(lm_shannon_vs_sheetwashfreq_male_log)
summary(anova_shannon_vs_sheetwashfreq_male_log)
TukeyHSD(anova_shannon_vs_sheetwashfreq_male_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_male_plot_shannon.png"
       , gg_shannon_sheetwashfreq_male
       , height=4, width=6)

#####Observed##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_observed_sheetwashfreq_male <- ggplot(samp_dat_wdiv_male, aes(x=`sheetwashfreq_binned`, y=Observed)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Observed Metrics for Male Sheet Wash Frequency")
gg_observed_sheetwashfreq_male

#ANOVA
# Set up our linear model 
lm_observed_vs_sheetwashfreq_male <- lm(Observed ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_male)
# Calculate AOV
anova_observed_vs_sheetwashfreq_male <- aov(lm_observed_vs_sheetwashfreq_male)
# Summarize to determine if there are significant differences
summary(anova_observed_vs_sheetwashfreq_male)
# Determine which groups are significant
tukey_sum_observed_male <- TukeyHSD(anova_observed_vs_sheetwashfreq_male)

# Kruskal-wallis test
kruskal_observed_male <- kruskal.test(Observed ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_male)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_observed_vs_sheetwashfreq_male_log <- lm(log(Observed) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_male)
anova_observed_vs_sheetwashfreq_male_log <- aov(lm_observed_vs_sheetwashfreq_male_log)
summary(anova_observed_vs_sheetwashfreq_male_log)
TukeyHSD(anova_observed_vs_sheetwashfreq_male_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_male_plot_observed.png"
       , gg_observed_sheetwashfreq_male
       , height=4, width=6)

#####Chao1##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_chao1_sheetwashfreq_male <- ggplot(samp_dat_wdiv_male, aes(x=`sheetwashfreq_binned`, y=Chao1)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Chao1 Metrics for Male Sheet Wash Frequency")
gg_chao1_sheetwashfreq_male

#ANOVA
# Set up our linear model 
lm_chao1_vs_sheetwashfreq_male <- lm(Chao1 ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_male)
# Calculate AOV
anova_chao1_vs_sheetwashfreq_male <- aov(lm_chao1_vs_sheetwashfreq_male)
# Summarize to determine if there are significant differences
summary(anova_chao1_vs_sheetwashfreq_male)
# Determine which groups are significant
tukey_sum_chao1_male <- TukeyHSD(anova_chao1_vs_sheetwashfreq_male)

# Kruskal-wallis test
kruskal_chao1_male <- kruskal.test(Chao1 ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_male)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_chao1_vs_sheetwashfreq_male_log <- lm(log(Chao1) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_male)
anova_chao1_vs_sheetwashfreq_male_log <- aov(lm_chao1_vs_sheetwashfreq_male_log)
summary(anova_chao1_vs_sheetwashfreq_male_log)
TukeyHSD(anova_chao1_vs_sheetwashfreq_male_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_male_plot_chao1.png"
       , gg_chao1_sheetwashfreq_male
       , height=4, width=6)

#####ACE##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_ace_sheetwashfreq_male <- ggplot(samp_dat_wdiv_male, aes(x=`sheetwashfreq_binned`, y=ACE)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "ACE Metrics for Male Sheet Wash Frequency")
gg_ace_sheetwashfreq_male

#ANOVA
# Set up our linear model 
lm_ace_vs_sheetwashfreq_male <- lm(ACE ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_male)
# Calculate AOV
anova_ace_vs_sheetwashfreq_male <- aov(lm_ace_vs_sheetwashfreq_male)
# Summarize to determine if there are significant differences
summary(anova_ace_vs_sheetwashfreq_male)
# Determine which groups are significant
tukey_sum_ace_male <- TukeyHSD(anova_ace_vs_sheetwashfreq_male)

# Kruskal-wallis test
kruskal_ace_male <- kruskal.test(ACE ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_male)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_ace_vs_sheetwashfreq_male_log <- lm(log(ACE) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_male)
anova_ace_vs_sheetwashfreq_male_log <- aov(lm_ace_vs_sheetwashfreq_male_log)
summary(anova_ace_vs_sheetwashfreq_male_log)
TukeyHSD(anova_ace_vs_sheetwashfreq_male_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_male_plot_ace.png"
       , gg_ace_sheetwashfreq_male
       , height=4, width=6)

#####Simpson##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_simpson_sheetwashfreq_male <- ggplot(samp_dat_wdiv_male, aes(x=`sheetwashfreq_binned`, y=Simpson)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Simpson Metrics for Male Sheet Wash Frequency")
gg_simpson_sheetwashfreq_male

#ANOVA
# Set up our linear model 
lm_simpson_vs_sheetwashfreq_male <- lm(Simpson ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_male)
# Calculate AOV
anova_simpson_vs_sheetwashfreq_male <- aov(lm_simpson_vs_sheetwashfreq_male)
# Summarize to determine if there are significant differences
summary(anova_simpson_vs_sheetwashfreq_male)
# Determine which groups are significant
tukey_sum_simpson_male <- TukeyHSD(anova_simpson_vs_sheetwashfreq_male)

# Kruskal-wallis test
kruskal_simpson_male <- kruskal.test(Simpson ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_male)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_simpson_vs_sheetwashfreq_male_log <- lm(log(Simpson) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_male)
anova_simpson_vs_sheetwashfreq_male_log <- aov(lm_simpson_vs_sheetwashfreq_male_log)
summary(anova_simpson_vs_sheetwashfreq_male_log)
TukeyHSD(anova_simpson_vs_sheetwashfreq_male_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_male_plot_simpson.png"
       , gg_simpson_sheetwashfreq_male
       , height=4, width=6)

#####InvSimpson##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_invsimpson_sheetwashfreq_male <- ggplot(samp_dat_wdiv_male, aes(x=`sheetwashfreq_binned`, y=InvSimpson)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "InvSimpson Metrics for Male Sheet Wash Frequency")
gg_invsimpson_sheetwashfreq_male

#ANOVA
# Set up our linear model 
lm_invsimp_vs_sheetwashfreq_male <- lm(InvSimpson ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_male)
# Calculate AOV
anova_invsimp_vs_sheetwashfreq_male <- aov(lm_invsimp_vs_sheetwashfreq_male)
# Summarize to determine if there are significant differences
summary(anova_invsimp_vs_sheetwashfreq_male)
# Determine which groups are significant
tukey_sum_invsimp_male <- TukeyHSD(anova_invsimp_vs_sheetwashfreq_male)

# Kruskal-wallis test
kruskal_invsimp_male <- kruskal.test(InvSimpson ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_male)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_invsimp_vs_sheetwashfreq_male_log <- lm(log(InvSimpson) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_male)
anova_invsimp_vs_sheetwashfreq_male_log <- aov(lm_invsimp_vs_sheetwashfreq_male_log)
summary(anova_invsimp_vs_sheetwashfreq_male_log)
TukeyHSD(anova_invsimp_vs_sheetwashfreq_male_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_male_plot_invsimpson.png"
       , gg_invsimpson_sheetwashfreq_male
       , height=4, width=6)

#####Fisher##### 
# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_fisher_sheetwashfreq_male <- ggplot(samp_dat_wdiv_male, aes(x=`sheetwashfreq_binned`, y=Fisher)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("low", "medium", "high")) +
  xlab("Sheet Wash Frequency") +
  labs(title = "Fisher Metrics for Male Sheet Wash Frequency")
gg_fisher_sheetwashfreq_male

#ANOVA
# Set up our linear model 
lm_fisher_vs_sheetwashfreq_male <- lm(Fisher ~ `sheetwashfreq_binned`, dat=samp_dat_wdiv_male)
# Calculate AOV
anova_fisher_vs_sheetwashfreq_male <- aov(lm_fisher_vs_sheetwashfreq_male)
# Summarize to determine if there are significant differences
summary(anova_fisher_vs_sheetwashfreq_male)
# Determine which groups are significant
tukey_sum_fisher_male <- TukeyHSD(anova_fisher_vs_sheetwashfreq_male)

# Kruskal-wallis test
kruskal_fisher_male <- kruskal.test(Fisher ~ `sheetwashfreq_binned`, data = samp_dat_wdiv_male)
# log the data and run ANOVA again to see that you also get better significance than just the ANOVA without the transformation
lm_fisher_vs_sheetwashfreq_male_log <- lm(log(Fisher) ~ `sheetwashfreq_binned`, data=samp_dat_wdiv_male)
anova_fisher_vs_sheetwashfreq_male_log <- aov(lm_fisher_vs_sheetwashfreq_male_log)
summary(anova_fisher_vs_sheetwashfreq_male_log)
TukeyHSD(anova_fisher_vs_sheetwashfreq_male_log)

#Results: no significance found

# save plot file 
ggsave(filename = "sheetwashfreq_male_plot_fisher.png"
       , gg_fisher_sheetwashfreq_male
       , height=4, width=6)

#####Phylogenetic diversity#####
# calculate Faith's phylogenetic diversity as PD
phylo_dist_sheetwashfreq_male <- pd(t(otu_table(dorms_rare_sheetwashfreq_male)), phy_tree(dorms_rare_sheetwashfreq_male),
                                    include.root=F)
# add PD to metadata table
sample_data(dorms_rare_sheetwashfreq_male)$PD <- phylo_dist_sheetwashfreq_male$PD
# plot any metadata category against the PD
plot.pd_sheetwashfreq_male <- ggplot(sample_data(dorms_rare_sheetwashfreq_male), aes(sheetwashfreq_binned, PD)) +
  geom_boxplot() +
  xlab("Sheet Washing Frequency") +
  ylab("Phylogenetic Diversity") +
  labs(title = "Phylogenetic Diversity for Male Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "medium", "high"))
# view plot
plot.pd_sheetwashfreq_male

# save plot file 
ggsave(filename = "sheetwashfreq_male_plot_pd.png"
       , plot.pd_sheetwashfreq_male
       , height=4, width=6)
