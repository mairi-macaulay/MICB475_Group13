####ALPHA DIVERSITY (1A - Sheet Washing Frequency)####

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

##Load in RData
load("dorms_rare.RData")
load("dorms_final.RData")

#remove NAs
dorms_rare <- subset_samples(dorms_rare, !is.na(sheetwashfreq_binned))

#### Alpha diversity ######
#view all metrics
plot_richness(dorms_rare)

# define which of the measures you want to show
plot_richness(dorms_rare, measures = c("Shannon","Chao1"))

# save into object, define x-axis as sheetwashfreq_binned and relabel, add a boxplot
gg_richness <- plot_richness(dorms_rare, x = "sheetwashfreq_binned", measures = c("Shannon","Chao1")) +
  xlab("Sheet Wash Frequency") +
  geom_boxplot()
gg_richness

# save plot file 
ggsave(filename = "sheetwashfreq_plot_richness.png"
       , gg_richness
       , height=4, width=6)

#tells you all the measures 
estimate_richness(dorms_rare)

## Phylogenetic diversity
# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(dorms_rare)), phy_tree(dorms_rare),
                 include.root=F)
# add PD to metadata table
sample_data(dorms_rare)$PD <- phylo_dist$PD
# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(dorms_rare), aes(sheetwashfreq_binned, PD)) +
  geom_boxplot() +
  xlab("Sheet Washing Frequency") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd