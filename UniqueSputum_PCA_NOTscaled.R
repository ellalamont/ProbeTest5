# PCA plot Sputum from ProbeTests 3,4,5

# 3/13/25 - making it with the NOT scaled


# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

# source("Import_data.R") # to get UniqueSputum_tpm and UniqueSputum_metadata
# Also UniqueSputum_tpm_NOTscaled

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


###########################################################
################ ALL UNIQUE SPUTUM PCA ####################
# Include all unique sputum samples, even those with low read counts

# Think I need to transform the data first
# UniqueSputum_tpm_t <- as.data.frame(t(UniqueSputum_tpm_NOTscaled %>% select(-Gene)))
UniqueSputum_tpm_t <- as.data.frame(t(UniqueSputum_tpm_NOTscaled))

# Remove columns that are all zero so the scale works for prcomp
UniqueSputum_tpm_t2 <- UniqueSputum_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(UniqueSputum_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 33.9% of variance
summary_PCA[2,1] # PC2 explains 22.5% of variance
summary_PCA[3,1] # PC3 explains 14.5% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, UniqueSputum_metadata, by = "SampleID")

fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Week, shape = Week)) + 
  # ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = Strain, text = Replicate)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8) +
  # geom_point(aes(fill = Week, shape = Week), size = 3, alpha = 0.8, stroke = 0.5) + # For thumbnail
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "NOT scaled PCA plot All Unique Sputum",
       subtitle = "All normal Depletion, no thresholds",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)

ggsave(fig_PC1vsPC2,
       file = "PCA_UniqueSputum_PC1vsPC2_NOTscaled.pdf",
       path = "PCA_Figures",
       width = 7, height = 5, units = "in")

# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Week# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")



###########################################################
############## >1M READS UNIQUE SPUTUM PCA ################
# Filter based on >1M read count (which is also at least 80% genes with at least 10 reads)

# Get the names of the samples with >1M reads
Unique_Sputum_1Mreads <- UniqueSputum_metadata %>% filter(N_Genomic > 1000000) %>% pull(SampleID)

# Filter to just have the samples >1M reads
UniqueSputum_1Mreads_tpm_t <- UniqueSputum_tpm_t %>% filter(row.names(UniqueSputum_tpm_t) %in% Unique_Sputum_1Mreads)

# Remove columns that are all zero so the scale works for prcomp
UniqueSputum_1Mreads_tpm_t2 <- UniqueSputum_1Mreads_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(UniqueSputum_1Mreads_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 36.0% of variance
summary_PCA[2,1] # PC2 explains 29.6% of variance
summary_PCA[3,1] # PC3 explains 22.5% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, UniqueSputum_metadata, by = "SampleID")

fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Week, shape = Week)) + 
  # ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = Strain, text = Replicate)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8) +
  # geom_point(aes(fill = Week, shape = Week), size = 3, alpha = 0.8, stroke = 0.5) + # For thumbnail
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  # geom_text_repel(aes(label = Run), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "NOT scaled PCA plot Unique Sputum",
       subtitle = "Sputum samples with >1M reads (Normal depletion)",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)

ggsave(fig_PC1vsPC2,
       file = "PCA_UniqueSputum_1Mreads_PC1vsPC2_v2.pdf",
       path = "PCA_Figures",
       width = 7, height = 5, units = "in")

# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Week# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")
