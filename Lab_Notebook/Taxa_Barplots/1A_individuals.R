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

# Joining OTU and metadata
otu_meta <- inner_join(metadata, otu_table, by = "ID")

# Create a directory to save the plot
dir.create("individual_plot", showWarnings = FALSE)

# Looping through each severity index to create a relative abundance measure
for (severity_level in unique(otu_meta$sheetwashfreq_binned)) {
  # Filter data for the current severity level
  df <- filter(otu_meta, sheetwashfreq_binned == severity_level)
  
  # Create a plot for the severity level
  p <- ggplot(data = df, aes(x = ID, y = abundance, fill = sheetwashfreq_binned)) +
    geom_col(color = "black") +
    labs(y = "Relative abundance", x = "Individuals", title = paste("Sheet Wash Frequency:", severity_level)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90))
  
  # Reset the graphics device
  dev.off()
  
  # Print the plot within RStudio
  print(p)
}