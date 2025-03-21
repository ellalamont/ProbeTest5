# PCA plot Sputum from ProbeTests 3,4,5 and the 3 captured broth samples from ProbeTest5
# 2/24/25


# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

# source("Import_data.R") # to get All_tpm

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )

poster_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=20), 
        axis.text.x = element_text(angle = 0, size=20, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )



###########################################################
################## ALL ProbeTest5 PCA #####################

ProbeTest5_tpm_2 <- ProbeTest5_tpm
# Adjust the names so they are slightly shorter
names(ProbeTest5_tpm_2) <- gsub(x = names(ProbeTest5_tpm_2), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# Grab the metadata I added to my_pipeSummary
ProbeTest5_metadata <- All_pipeSummary %>% select(1, 13:24) %>%
  filter(Run == "ProbeTest5")

# add rownames to the tpm and metadata dataframes
rownames(ProbeTest5_tpm_2) <- ProbeTest5_tpm_2[,1] # add the rownames
ProbeTest5_tpm_2 <- ProbeTest5_tpm_2[,-1] # Remove the old column of rownames
# rownames(my_metadata) <- my_metadata[,1] # add the rownames
# my_metadata <- my_metadata[,-1] # Remove the old column of rownames

# Think I need to transform the data first
ProbeTest5_tpm_2_t <- as.data.frame(t(ProbeTest5_tpm_2))

# Remove columns that are all zero so the scale works for prcomp
ProbeTest5_tpm_2_t2 <- ProbeTest5_tpm_2_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(ProbeTest5_tpm_2_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 31.3% of variance
summary_PCA[2,1] # PC2 explains 7.1% of variance
summary_PCA[3,1] # PC3 explains 6.6% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, ProbeTest5_metadata, by = "SampleID", )

fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = Probe)) + 
  # ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = Strain, text = Replicate)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8) +
  # geom_point(aes(fill = Week, shape = Week), size = 3, alpha = 0.8, stroke = 0.5) + # For thumbnail
  # scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  # scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA plot All ProbeTest5",
       subtitle = "All normal Depletion, no thresholds",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Week# , 
                  # colors = c12,
                  # text = ~Replicate
)
# PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")


###########################################################
################### SPUTUM WITH BROTH #####################

# Start with All_tpm
# Select all the broth and all the sputum that are >1M reads
mySubset_tpm <- All_tpm %>% select(contains("Broth") | all_of(Unique_Sputum_1Mreads))

# Transform the data
mySubset_tpm_t <- as.data.frame(t(mySubset_tpm))

# Remove columns that are all zero so the scale works for prcomp
mySubset_tpm_t2 <- mySubset_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(mySubset_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 33.8% of variance
summary_PCA[2,1] # PC2 explains 23.2% of variance
summary_PCA[3,1] # PC3 explains 15.3% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_metadata, by = "SampleID", )

fig_PC1vsPC2 <- my_PCA_df %>% 
  mutate(Labelling = c("Captured broth", "Captured broth", "Captured broth", "Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W0 sputum", "W0 sputum", "W2 sputum", "W2 sputum", "W2 sputum")) %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Labelling, shape = Labelling), size = 6, alpha = 0.8, stroke = 0.8) + 
  scale_fill_manual(values=c(`W0 sputum` = "#0072B2", `W2 sputum` = "#E66900", `Captured broth`= "maroon", `Not captured broth`= "#999999")) +  
  scale_shape_manual(values=c(`W0 sputum` = 21, `W2 sputum` = 22, `Captured broth`= 23, `Not captured broth`= 23)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA Unique sputum >1M reads with captured and not captured broth",
       # subtitle = "All normal Depletion, no thresholds",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)
ggsave(fig_PC1vsPC2,
       file = "PCA_UniqueSputum1Mreads_With_Broth.pdf",
       path = "PCA_Figures",
       width = 8, height = 5, units = "in")


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- my_PCA_df %>% 
  mutate(Labelling = c("Captured broth", "Captured broth", "Captured broth", "Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W0 sputum", "W0 sputum", "W2 sputum", "W2 sputum", "W2 sputum")) %>%
  plot_ly(x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Labelling# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")



###########################################################
################ SPIKE THP1 WITH BROTH ####################

# Compare all the THP1 spiked samples with the broth samples (captured and not). These aren't the same culture but should be similar. Which broth is truth??

# Need to subset to only keep the samples with >1M reads, so will lose some of the 1e2 and 1e3 samples
Samples_1Mread_list <- All_pipeSummary %>% 
  filter(N_Genomic > 1000000) %>% 
  pull(SampleID)

# Still keeping samples that dont have 80% genes with at least 10 reads mapping.. Just keep those...
Samples_1Mread.80PGenes_list <- All_pipeSummary %>% 
  filter(N_Genomic > 1000000) %>% 
  filter(AtLeast.10.Reads >= 4499*0.8) %>% 
  pull(SampleID)

# Start with ProbeTest5_tpm because I only want the ones in this run
# Select all the broth and all the sputum that are >1M reads
mySubset_tpm <- ProbeTest5_tpm_2 %>% 
  select(any_of(Samples_1Mread.80PGenes_list)) %>% # only keep samples with >1M reads
  select(!contains("S_")) # remove the sputum 
  

# Think I need to transform the data first
mySubset_tpm_t <- as.data.frame(t(mySubset_tpm))

# Remove columns that are all zero so the scale works for prcomp
mySubset_tpm_t2 <- mySubset_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(mySubset_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 62.9% of variance
summary_PCA[2,1] # PC2 explains 17.4% of variance
summary_PCA[3,1] # PC3 explains 8.0% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, ProbeTest5_metadata, by = "SampleID", )

fig_PC1vsPC2 <- my_PCA_df %>%
  # mutate(Labelling = c("Captured broth", "Captured broth", "Captured broth", "Not captured broth", "Not captured broth", "Not captured broth", "1e3", "1e3", "1e4", "1e4", "1e4", "1e5", "1e5", "1e5", "1e6", "1e6", "1e6", "1e6")) %>% # This for just using the 1M threshold, not the 80% genes threshold
  mutate(Labelling = c("Captured broth", "Captured broth", "Captured broth", "Not captured broth", "Not captured broth", "Not captured broth", "1e5", "1e5", "1e5", "1e6", "1e6", "1e6", "1e6")) %>%
  ggplot(aes(x = PC1, y = PC2, text = SampleID)) + 
  geom_point(aes(fill = Labelling, shape = Labelling), size = 6, alpha = 0.8, stroke = 0.8) + 
  scale_fill_manual(values=c(`Captured broth`= "maroon", `Not captured broth`= "#999999", `1e3` = "#B2EBF2", `1e4` = "#81D4FA", `1e5` = "#0288D1", `1e6` = "#01579B")) +  
  scale_shape_manual(values=c(`Captured broth`= 23, `Not captured broth`= 23, `1e3` = 21, `1e4` = 21, `1e5` = 21, `1e6` = 21)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA plot All ProbeTest5",
       subtitle = "All >1 Million reads and >=80% genes with at least 10 reads aligning ",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)

ggsave(fig_PC1vsPC2,
       file = "PCA_THP1Spiked1MReads.80PGenes_With_Broth.pdf",
       path = "PCA_Figures",
       width = 8, height = 5, units = "in")


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- my_PCA_df %>%
  # mutate(Labelling = c("Captured broth", "Captured broth", "Captured broth", "Not captured broth", "Not captured broth", "Not captured broth", "1e3", "1e3", "1e4", "1e4", "1e4", "1e5", "1e5", "1e5", "1e6", "1e6", "1e6", "1e6")) %>%
  mutate(Labelling = c("Captured broth", "Captured broth", "Captured broth", "Not captured broth", "Not captured broth", "Not captured broth", "1e5", "1e5", "1e5", "1e6", "1e6", "1e6", "1e6")) %>%
  plot_ly(x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Labelling# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")

###########################################################
############### NOT scaled SPUTUM WITH BROTH ##############

# 3/13/25 : Switching this to be NOT scaled TPM and uncaptured broth only

# Start with All_tpm_NOTscaled
# Select all the broth and all the sputum that are >1M reads
mySubset_tpm_NOTscaled <- All_tpm_NOTscaled %>% select(c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6") | all_of(Unique_Sputum_1Mreads))

# Transform the data
mySubset_tpm_t <- as.data.frame(t(mySubset_tpm_NOTscaled))

# Remove columns that are all zero so the scale works for prcomp
mySubset_tpm_t2 <- mySubset_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(mySubset_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 32.9% of variance
summary_PCA[2,1] # PC2 explains 24.1% of variance
summary_PCA[3,1] # PC3 explains 18.7% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_metadata, by = "SampleID", )

fig_PC1vsPC2 <- my_PCA_df %>% 
  # mutate(Labelling = c("Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W2 sputum", "W2 sputum", "W0 sputum", "W0 sputum", "W2 sputum")) %>% # NEED TO CHECK BY HAND THAT THESE ARE CORRECT LABELS!
  mutate(Labelling = c("Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W0 sputum", "W0 sputum", "W2 sputum", "W2 sputum", "W2 sputum")) %>% # NEED TO CHECK BY HAND THAT THESE ARE CORRECT LABELS!
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Labelling, shape = Labelling), size = 6, alpha = 0.8, stroke = 0.8) + 
  scale_fill_manual(values=c(`W0 sputum` = "#0072B2", `W2 sputum` = "#E66900", `Not captured broth`= "#999999")) +  
  scale_shape_manual(values=c(`W0 sputum` = 21, `W2 sputum` = 22, `Not captured broth`= 23)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA Unique sputum >1M reads with captured and not captured broth",
       # subtitle = "All normal Depletion, no thresholds",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)
ggsave(fig_PC1vsPC2,
       file = "PCA_UniqueSputum1Mreads_With_Broth_v2.pdf",
       path = "PCA_Figures",
       width = 8, height = 5, units = "in")


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- my_PCA_df %>% 
  mutate(Labelling = c("Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W2 sputum", "W2 sputum", "W0 sputum", "W0 sputum", "W2 sputum")) %>%
  plot_ly(x = ~PC1, y = ~PC2, z = ~PC3,
          type = "scatter3d", mode = "markers",
          color = ~Labelling# , 
          # colors = c12,
          # text = ~Replicate
  )
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")



# For poster
fig_PC1vsPC2_poster <- my_PCA_df %>% 
  # mutate(Labelling = c("Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W2 sputum", "W2 sputum", "W0 sputum", "W0 sputum", "W2 sputum")) %>% # NEED TO CHECK BY HAND THAT THESE ARE CORRECT LABELS!
  mutate(Labelling = c("Bacterial broth", "Bacterial broth", "Bacterial broth", "Week 0 sputum", "Week 0 sputum", "Week 0 sputum", "Week 2 sputum", "Week 2 sputum", "Week 2 sputum")) %>% # NEED TO CHECK BY HAND THAT THESE ARE CORRECT LABELS!
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Labelling, shape = Labelling), size = 6, alpha = 0.8, stroke = 0.8) + 
  scale_fill_manual(values=c(`Week 0 sputum` = "#0072B2", `Week 2 sputum` = "#E66900", `Bacterial broth`= "#999999")) +  
  scale_shape_manual(values=c(`Week 0 sputum` = 21, `Week 2 sputum` = 22, `Bacterial broth`= 23)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = NULL,
       # subtitle = "All normal Depletion, no thresholds",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  poster_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2_poster
# ggplotly(fig_PC1vsPC2)
ggsave(fig_PC1vsPC2_poster,
       file = "PCA_UniqueSputum1Mreads_With_Broth_v1.pdf",
       path = "Poster_Figures",
       width = 7, height = 4.5, units = "in")







###########################################################
##################### W0 AND BROTH PCA ####################

# 3/20/25 Post Group meeting: What does the PCA look like without W2

# Start with All_tpm_NOTscaled
# Select all the broth and the W0 sputum samples 
# W0 samples: "S_250754", "S_355466", "S_503557" 
mySubset_tpm_NOTscaled <- All_tpm_NOTscaled %>% select(c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6", "S_250754", "S_355466", "S_503557" ))

# Transform the data
mySubset_tpm_t <- as.data.frame(t(mySubset_tpm_NOTscaled))

# Remove columns that are all zero so the scale works for prcomp
mySubset_tpm_t2 <- mySubset_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(mySubset_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 50.1% of variance
summary_PCA[2,1] # PC2 explains 23.3% of variance
summary_PCA[3,1] # PC3 explains 19.5% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_metadata, by = "SampleID", )

fig_PC1vsPC2 <- my_PCA_df %>% 
  mutate(Labelling = c("Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W0 sputum", "W0 sputum")) %>% # NEED TO CHECK BY HAND THAT THESE ARE CORRECT LABELS!
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Labelling, shape = Labelling), size = 6, alpha = 0.8, stroke = 0.8) + 
  scale_fill_manual(values=c(`W0 sputum` = "#0072B2", `Not captured broth`= "#999999")) +  
  scale_shape_manual(values=c(`W0 sputum` = 21, `W2 sputum` = 22, `Not captured broth`= 23)) + 
  # geom_text_repel(aes(label = Run), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA W0 sputum >1M reads with not captured broth",
       # subtitle = "All normal Depletion, no thresholds",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)
ggsave(fig_PC1vsPC2,
       file = "PCA_W0_With_Broth_v1.pdf",
       path = "PCA_Figures",
       width = 8, height = 5, units = "in")

# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- my_PCA_df %>% 
  mutate(Labelling = c("Not captured broth", "Not captured broth", "Not captured broth", "W0 sputum", "W0 sputum", "W0 sputum")) %>%
  plot_ly(x = ~PC1, y = ~PC2, z = ~PC3,
          type = "scatter3d", mode = "markers",
          color = ~Labelling# , 
          # colors = c12,
          # text = ~Replicate
  )
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")
