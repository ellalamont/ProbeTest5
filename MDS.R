# MDS plot with my W0, W2, and broth
# E. Lamont
# 5/1/25

# https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/122-multidimensional-scaling-essentials-algorithms-and-r-code/
# https://www.geeksforgeeks.org/multidimensional-scaling-using-r/



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
###################### PROCESS DATA #######################

mySubset_tpm_NOTscaled <- All_tpm_NOTscaled %>% select(c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6") | all_of(Unique_Sputum_1Mreads))

# Transform the data
mySubset_tpm_t <- as.data.frame(t(mySubset_tpm_NOTscaled))

# Remove columns that are all zero so the scale works for prcomp
mySubset_tpm_t2 <- mySubset_tpm_t %>% select_if(colSums(.) != 0)

###########################################################
##################### CLASSICAL MDS #######################

# Cmpute MDS
mds <- cmdscale(dist(mySubset_tpm_t2))
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
       subtitle = "Using TPM",
       x = "Dimension 1",
       y = "Dimension 2") +
  my_plot_themes
# my_plot_themes_thumbnail
fig_MDS
ggsave(fig_MDS,
       file = "MDS_UniqueSputum1Mreads_With_Broth_v1.pdf",
       path = "MDS_Figures",
       width = 8, height = 5, units = "in")


