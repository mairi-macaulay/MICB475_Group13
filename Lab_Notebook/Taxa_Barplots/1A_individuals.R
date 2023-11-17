library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)

# Load Rdata
load("../Phyloseq/dorms_final_sheetwashfreq.RData")

# Extracting OTU data
otu_table <- data.frame(t(otu_table(dorms_final)))
otu_table$ID <- rownames(otu_table)

# Extracting metadata
metadata <- data.frame(sample_data(dorms_final))
metadata$ID <- rownames(metadata)

# Load the raw taxonomy file
tax <- data.frame(tax_table(dorms_final))

# Formatting the taxa dataframe and cleaning names
tax_mat <- tax[,-1]
tax_mat <- data.frame(tax_mat)
tax_mat$Phylum <- gsub("^...","",tax_mat$Phylum)
tax_mat$Class <- gsub("^...","",tax_mat$Class)
tax_mat$Order <- gsub("^...","",tax_mat$Order)
tax_mat$Family <- gsub("^...","",tax_mat$Family)
tax_mat$Genus <- gsub("^...","",tax_mat$Genus)
tax_mat$Species <- gsub("^...","",tax_mat$Species)
tax_mat$ASV <- rownames(tax_mat)

# Joining OTU and metadata and taxanomic information
otu_meta <- inner_join(metadata, otu_table, by = "ID")

#transforming the OTU matrix to a single column called abundance. 53 represents the number of metadata columns we wnat to exlcude
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:53))

#Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")

#Generate graphs for 'low', 'medium', and 'high' sheet wash frequency
sheetwash_levels <- c("low", "medium", "high")

#list to store plots
plots <- list()

#Generate plots for each sheet wash freq level
for (level in sheetwash_levels) {
  #filter data for current sheet wash fre level
  current_data <- grouped_taxa %>%
    filter(sheetwashfreq_binned == level)
  
  #Plotting at phylum level
  plot <- ggplot(data = current_data, aes(x = ID, fill = Phylum)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = -90)) +
    labs(y = "Relative Abundance", x = "Sheet Wash Individuals") +
    theme_classic() +
    theme(axis.text = element_text(size = 5, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"))
  
  #Store plot in list
  plots[[level]] <- plot
}

#show plots
plots[['low']]
plots[['medium']]
plots[['high']]
