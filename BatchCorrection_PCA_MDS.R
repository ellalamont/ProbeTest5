# Batch Correction and MDS plot, following what Mark has done....
# E. Lamont
# 5/1/25

# Just looking at the Samples I have been moving forward with (>1M and uncaptured broth)

# BiocManager::install("sva")
library(sva) # For ComBat_seq batch correction
library(edgeR)

source("Import_data.R") # SputumBroth_RawReads

###########################################################
###################### PROCESS DATA #######################

# Mark starts with READS_M for this
# SputumBroth_RawReads

SputumBroth_RawReads$Gene <- NULL
SputumBroth_RawReads <- as.matrix(SputumBroth_RawReads)

# Subset the metadata
SputumBroth_pipeSummary <- All_pipeSummary %>% filter(SampleID %in% c(Unique_Sputum_1Mreads, "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")) %>%
  mutate(Batch = case_when(
    Run == "ProbeTest4" ~ 1, 
    Run == "ProbeTest5" ~ 2
  ))

###########################################################
#################### BATCH CORRECTION #####################
# This basically all taken from Mark and I'm not sure what's happening

# https://academic.oup.com/nargab/article/2/3/lqaa078/5909519?login=true

batch <- SputumBroth_pipeSummary$Batch # This isn't in the same order as the columns so I don't think this is right...
# Just do it manually
batch <- c(2,2,2,1,1,1,2,2,2)
counts_corrected <- ComBat_seq(SputumBroth_RawReads, batch = batch)


###########################################################
########################### PCA ###########################
# Just do the PCA the way I have been doing it

# Convert the batch corrected counts to counts per million (cpm)
merged_cpm <- cpm(counts_corrected)
# merged_cpm_log2 <- cpm(counts_corrected, log = T)

# transform the data 
merged_cpm_t <- as.data.frame(t(merged_cpm))
# merged_cpm_log2_t <- as.data.frame(t(merged_cpm_log2))

# Remove columns that are all zero so the scale works for prcomp
merged_cpm_t <- merged_cpm_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 
# merged_cpm_log2_t <- merged_cpm_log2_t %>% select_if(colSums(.) != 0) # Stays at 4499 genes

# Make the actual PCA
# Need to choose here if using log2 or not
my_PCA <- prcomp(merged_cpm_t, scale = TRUE)
# my_PCA <- prcomp(merged_cpm_log2_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 34.4% of variance
summary_PCA[2,1] # PC2 explains 25.2% of variance
summary_PCA[3,1] # PC3 explains 18.4% of variance

# Add the metadata
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_metadata, by = "SampleID") %>% 
  mutate(Labelling = c("Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W0 sputum", "W0 sputum", "W2 sputum", "W2 sputum", "W2 sputum")) # NEED TO CHECK BY HAND THAT THESE ARE CORRECT LABELS!

# MAKE PCA PLOT with GGPLOT 
fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Labelling, shape = Labelling)) + 
  geom_point(size = 6, alpha = 0.9, stroke = 0.8) +
  scale_fill_manual(values=c(`W0 sputum` = "#0072B2", `W2 sputum` = "#E66900", `Not captured broth`= "#999999")) +  
  scale_shape_manual(values=c(`W0 sputum` = 21, `W2 sputum` = 22, `Not captured broth`= 23)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: Unique sputum >1M reads with not captured broth",
       subtitle = "Batch corrected raw reads -> CPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2
ggsave(fig_PC1vsPC2,
       file = "PCA_UniqueSputum1Mreads_BatchCorrect.CPM_v1.pdf",
       path = "PCA_Figures",
       width = 8, height = 5, units = "in")


###########################################################
##################### CLASSICAL MDS #######################

# Convert the batch corrected counts to counts per million (cpm)
merged_cpm <- cpm(counts_corrected)
# merged_cpm_log2 <- cpm(counts_corrected, log = T)

# transform the data 
merged_cpm_t <- as.data.frame(t(merged_cpm))
# merged_cpm_log2_t <- as.data.frame(t(merged_cpm_log2))

# Remove columns that are all zero so the scale works for prcomp
merged_cpm_t <- merged_cpm_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 
# merged_cpm_log2_t <- merged_cpm_log2_t %>% select_if(colSums(.) != 0) # Stays at 4499 genes

# Cmpute MDS
mds <- cmdscale(dist(merged_cpm_t))
colnames(mds) <- c("Dim.1", "Dim.2")

# Add the metadata
my_MDS_df <- as.data.frame(mds) %>%
  rownames_to_column(var = "SampleID")
my_MDS_df <- merge(my_MDS_df, my_metadata, by = "SampleID") %>%
  mutate(Labelling = c("Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W0 sputum", "W0 sputum", "W2 sputum", "W2 sputum", "W2 sputum")) # NEED TO CHECK BY HAND THAT THESE ARE CORRECT LABELS!


fig_MDS <- my_MDS_df %>%
  ggplot(aes(x = Dim.1, y = Dim.2, fill = Labelling, shape = Labelling)) + 
  geom_point(size = 6, alpha = 0.9, stroke = 0.8) +
  scale_fill_manual(values=c(`W0 sputum` = "#0072B2", `W2 sputum` = "#E66900", `Not captured broth`= "#999999")) +  
  scale_shape_manual(values=c(`W0 sputum` = 21, `W2 sputum` = 22, `Not captured broth`= 23)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "MDS: Unique sputum >1M reads with not captured broth",
       subtitle = "Batch corrected raw reads -> CPM",
       x = "Dimension 1",
       y = "Dimension 2") +
  my_plot_themes
# my_plot_themes_thumbnail
fig_MDS
ggsave(fig_MDS,
       file = "MDS_UniqueSputum1Mreads_BatchCorrect.CPM_v1.pdf",
       path = "MDS_Figures",
       width = 8, height = 5, units = "in")







