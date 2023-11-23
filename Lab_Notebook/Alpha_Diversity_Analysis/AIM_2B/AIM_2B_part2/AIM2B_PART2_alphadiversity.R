####ALPHA DIVERSITY (2B - Sheet Washing Frequency and Gender)####

#set seed
set.seed(1)

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)

####COMPARE FEMALE AND MALE - Two-Way ANOVAs####

##Load in RData
load("../../../Phyloseq/dorms_rare_sheetwashfreq.RData")

samp_dat_wdiv <- data.frame(sample_data(dorms_rare), estimate_richness(dorms_rare))

##Remove medium sheet washing frequency
samp_dat_wdiv_filt <- samp_dat_wdiv[samp_dat_wdiv$sheetwashfreq_binned != "medium", ]

#####Shannon#####
# plot the sheetwashfreq and gender against Shannon
plot_shannon_gender_sheetwashfreq <- ggplot(samp_dat_wdiv_filt) + geom_boxplot(aes(x=sheetwashfreq_binned, y=Shannon)) +
  facet_grid(~factor(`sex`, levels=c("female","male"))) +
  xlab("Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "high"))
plot_shannon_gender_sheetwashfreq

# run the 2-way ANOVA
shannon_sheetwashfreq_gender <- lm(Shannon ~ `sex`*`sheetwashfreq_binned`, data=samp_dat_wdiv_filt)
summary(aov(shannon_sheetwashfreq_gender))
TukeyHSD(aov(shannon_sheetwashfreq_gender))

#Results: no significance found

#save plot
ggsave(filename = "2B_part2_plot_shannon.png"
       , plot_shannon_gender_sheetwashfreq
       , height=4, width=6)

#####Observed#####
# plot the sheetwashfreq and gender against Observed
plot_observed_gender_sheetwashfreq <- ggplot(samp_dat_wdiv_filt) + geom_boxplot(aes(x=sheetwashfreq_binned, y=Observed)) +
  facet_grid(~factor(`sex`, levels=c("female","male"))) +
  xlab("Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "high"))
plot_observed_gender_sheetwashfreq

# run the 2-way ANOVA
observed_sheetwashfreq_gender <- lm(Observed ~ `sex`*`sheetwashfreq_binned`, data=samp_dat_wdiv_filt)
summary(aov(observed_sheetwashfreq_gender))
TukeyHSD(aov(observed_sheetwashfreq_gender))

#Results: no significance found

#save plot
ggsave(filename = "2B_part2_plot_observed.png"
       , plot_observed_gender_sheetwashfreq
       , height=4, width=6)

#####Chao1#####
# plot the sheetwashfreq and gender against Chao1
plot_chao1_gender_sheetwashfreq <- ggplot(samp_dat_wdiv_filt) + geom_boxplot(aes(x=sheetwashfreq_binned, y=Chao1)) +
  facet_grid(~factor(`sex`, levels=c("female","male"))) +
  xlab("Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "high"))
plot_chao1_gender_sheetwashfreq

# run the 2-way ANOVA
chao1_sheetwashfreq_gender <- lm(Chao1 ~ `sex`*`sheetwashfreq_binned`, data=samp_dat_wdiv_filt)
summary(aov(chao1_sheetwashfreq_gender))
TukeyHSD(aov(chao1_sheetwashfreq_gender))

#Results: no significance found

#save plot
ggsave(filename = "2B_part2_plot_chao1.png"
       , plot_chao1_gender_sheetwashfreq
       , height=4, width=6)

#####ACE#####
# plot the sheetwashfreq and gender against ACE
plot_ace_gender_sheetwashfreq <- ggplot(samp_dat_wdiv_filt) + geom_boxplot(aes(x=sheetwashfreq_binned, y=ACE)) +
  facet_grid(~factor(`sex`, levels=c("female","male"))) +
  xlab("Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "high"))
plot_ace_gender_sheetwashfreq

# run the 2-way ANOVA
ace_sheetwashfreq_gender <- lm(ACE ~ `sex`*`sheetwashfreq_binned`, data=samp_dat_wdiv_filt)
summary(aov(ace_sheetwashfreq_gender))
TukeyHSD(aov(ace_sheetwashfreq_gender))

#Results: no significance found

#save plot
ggsave(filename = "2B_part2_plot_ace.png"
       , plot_ace_gender_sheetwashfreq
       , height=4, width=6)

#####Simpson#####
# plot the sheetwashfreq and gender against Simpson
plot_simpson_gender_sheetwashfreq <- ggplot(samp_dat_wdiv_filt) + geom_boxplot(aes(x=sheetwashfreq_binned, y=Simpson)) +
  facet_grid(~factor(`sex`, levels=c("female","male"))) +
  xlab("Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "high"))
plot_simpson_gender_sheetwashfreq

# run the 2-way ANOVA
simpson_sheetwashfreq_gender <- lm(Simpson ~ `sex`*`sheetwashfreq_binned`, data=samp_dat_wdiv_filt)
summary(aov(simpson_sheetwashfreq_gender))
TukeyHSD(aov(simpson_sheetwashfreq_gender))

#Results: no significance found

#save plot
ggsave(filename = "2B_part2_plot_simpson.png"
       , plot_simpson_gender_sheetwashfreq
       , height=4, width=6)

#####InvSimpson#####
# plot the sheetwashfreq and gender against InvSimpson
plot_invsimp_gender_sheetwashfreq <- ggplot(samp_dat_wdiv_filt) + geom_boxplot(aes(x=sheetwashfreq_binned, y=InvSimpson)) +
  facet_grid(~factor(`sex`, levels=c("female","male"))) +
  xlab("Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "high"))
plot_invsimp_gender_sheetwashfreq

# run the 2-way ANOVA
invsimp_sheetwashfreq_gender <- lm(InvSimpson ~ `sex`*`sheetwashfreq_binned`, data=samp_dat_wdiv_filt)
summary(aov(invsimp_sheetwashfreq_gender))
TukeyHSD(aov(invsimp_sheetwashfreq_gender))

#Results: no significance found

#save plot
ggsave(filename = "2B_part2_plot_invsimp.png"
       , plot_invsimp_gender_sheetwashfreq
       , height=4, width=6)

#####Fisher#####
# plot the sheetwashfreq and gender against Fisher
plot_fisher_gender_sheetwashfreq <- ggplot(samp_dat_wdiv_filt) + geom_boxplot(aes(x=sheetwashfreq_binned, y=Fisher)) +
  facet_grid(~factor(`sex`, levels=c("female","male"))) +
  xlab("Sheet Wash Frequency") +
  scale_x_discrete(limits = c("low", "high"))
plot_fisher_gender_sheetwashfreq

# run the 2-way ANOVA
fisher_sheetwashfreq_gender <- lm(Fisher ~ `sex`*`sheetwashfreq_binned`, data=samp_dat_wdiv_filt)
summary(aov(fisher_sheetwashfreq_gender))
TukeyHSD(aov(fisher_sheetwashfreq_gender))

#Results: no sigificance found

#save plot
ggsave(filename = "2B_part2_plot_fisher.png"
       , plot_fisher_gender_sheetwashfreq
       , height=4, width=6)