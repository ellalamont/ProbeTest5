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
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")


###########################################################
################# THP1 SPIKE WITH BROTH ###################

# Compare all the THP1 spiked samples with the broth samples (captured and not). These aren't the same culture but should be similar. Which broth is truth??

# add rownames to the tpm and metadata dataframes
rownames(ProbeTest5_tpm) <- ProbeTest5_tpm[,1] # add the rownames
ProbeTest5_tpm <- ProbeTest5_tpm[,-1] # Remove the old column of rownames
names(ProbeTest5_tpm) <- gsub(x = names(ProbeTest5_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

mySubset_tpm <- ProbeTest5_tpm %>% select(contains("Broth") | contains("S_"))

# Think I need to transform the data first
mySubset_tpm_t <- as.data.frame(t(mySubset_tpm))

# Remove columns that are all zero so the scale works for prcomp
mySubset_tpm_t2 <- mySubset_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(mySubset_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 47.5% of variance
summary_PCA[2,1] # PC2 explains 21.9% of variance
summary_PCA[3,1] # PC3 explains 18.7% of variance

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
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")

