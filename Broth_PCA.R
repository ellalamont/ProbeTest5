# PCA plot Sputum from ProbeTests 5
# 2/16/25

# source("Import_data.R") # to get Broth_tpm_CorrectScales and Broth_metadata

# Broth_tpm_CorrectScales
# This has the captured samples (1-3) that are scaled and the not captured samples (4-6) that are not scaled!

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above


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
#################### ALL BROTH PCA ########################

# Remove the gene column
tmp <- Broth_tpm_CorrectScales %>% select(-Gene)

# Think I need to transform the data first
Broth_tpm_t <- as.data.frame(t(tmp))

# Remove columns that are all zero so the scale works for prcomp
Broth_tpm_t2 <- Broth_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(Broth_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 89.5% of variance
summary_PCA[2,1] # PC2 explains 4.7% of variance
summary_PCA[3,1] # PC3 explains 2.6% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, Broth_metadata, by = "SampleID")

fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Probe, shape = Probe)) + 
  # ggplot(aes(x = PC1, y = PC2, color = Sample_Type, shape = Strain, text = Replicate)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8) +
  # geom_point(aes(fill = Week, shape = Week), size = 3, alpha = 0.8, stroke = 0.5) + # For thumbnail
  scale_fill_manual(values=c(`Captured` = "#E69F00", `Not captured` = "#999999")) +  
  scale_shape_manual(values=c(`Captured` = 21, `Not captured` = 22)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA plot H37Ra broth",
       subtitle = "Captured Scaled, Not captured Not scaled",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)

ggsave(fig_PC1vsPC2,
       file = "PCA_Broth_PC1vsPC2.pdf",
       path = "PCA_Figures/Broth",
       width = 7, height = 5, units = "in")

# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Probe# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")
