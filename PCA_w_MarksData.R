# PCA plot with Marks data, sputum, broth (uncaptured) etc
# 3/13/15

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
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
################### IMPORT MARK's DATA ####################

# Not sure if Mark's are scaled or not....

# Import the tpm
mark_tpm <- read.csv("DataFromMark/Data/Updated.Mtb.Expression.Gene.Data.TPM.csv")
# add rownames to the tpm and metadata dataframes
rownames(mark_tpm) <- mark_tpm[,1] # add the rownames
mark_tpm <- mark_tpm %>% rename(Gene = X)

# Import the metatdata
mark_metadata <- read.csv("DataFromMark/Data/Annotation_including_new_timecourse_batch.csv")

# I think I just need the GroupStrain the Mimic_timepoint and Batch
# GroupStrain should be the same as Sample_Type
colnames(mark_metadata)[colnames(mark_metadata) == "GroupStrain"] ="Sample_Type"
mark_metadata_2 <- mark_metadata %>% select(SampleID, Sample_Type, Batch, Mimic_timepoint)

# Just keep the mimic and rabbit samples
sample_list <- mark_metadata_2 %>%
  filter(Sample_Type %in% c("mimic", "rabbit")) %>%
  pull(SampleID)
mark_tpm2 <- mark_tpm %>% select(all_of(c("Gene", sample_list)))
mark_metadata3 <- mark_metadata_2 %>% filter(SampleID %in% sample_list)

###########################################################
#################### COMBINE THE TPMs #####################

# Get all of mine that are NOT scaled
# ProbeTest5_tpm_NOTscaled
# All_tpm_NOTscaled

# Get just my samples that I want to add
# All Unique sputum, uncaptured broth, captured THP1 spikes (1e6)
# Lets pull everything I need from ProbeTest5 and then grab the unique sputum

ProbeTest5_samples <- c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6", "THP1_1e6_1a", "THP1_1e6_2b", "THP1_1e6_3a", "Gene")

subset1 <- ProbeTest5_tpm_NOTscaled %>% select(all_of(ProbeTest5_samples))

Unique_Sputum <- c("S_250754", "S_355466", "S_503557", "S_349941_Probe_3D_25", "S_503937", "S_575533_MtbrRNA", "S_577208", "S_351946_Probe_4A_100", "S_687338_Probe_4A_100", "Gene")

subset2 <- All_tpm_NOTscaled %>% select(all_of(Unique_Sputum))

# Merge the TPMs
combined_TPM <- merge(subset1, subset2, by = "Gene", all = T)
combined_TPM <- merge(combined_TPM, mark_tpm2, by = "Gene", all = T)
rownames(combined_TPM) <- combined_TPM[,1] # add the rownames

# Merge the metadatas
my_metadata_subset <- my_metadata %>%
  filter(SampleID %in% c(ProbeTest5_samples, Unique_Sputum))
combined_metadata <- merge(my_metadata_subset, mark_metadata3, all = T)

# Merge Mimic_timepoint and Week
combined_metadata <- combined_metadata %>%
  mutate(Merged_Timpoint = coalesce(Mimic_timepoint, Week)) %>%
  mutate(Merged_Timpoint = na_if(Merged_Timpoint, "na"))


###########################################################
######################## MAKE PCA #########################

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

# Think I need to transform the data first
my_tpm_t <- as.data.frame(t(combined_TPM %>% select(-Gene)))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # Adding round and digits here to set number of digits after the decimal place.
summary_PCA[1,1] # PC1 explains 24.4% of variance
summary_PCA[2,1] # PC2 explains 15.7% of variance
summary_PCA[3,1] # PC3 explains 5.9% of variance


###########################################################
################ MAKE PCA PLOT with GGPLOT ################

my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, combined_metadata, by = "SampleID")

fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, color = Sample_Type, label = SampleID, label2 = Week, label3 = Batch)) + 
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Merged_Timpoint), color = "black", size = 2, max.overlaps = Inf) + 
  scale_color_manual(values = c(`Sputum` = "#0072B2", `THP1` = "#FF7F00", `rabbit` = "orchid1", `Broth` = "black", `mimic` = "purple")) + 
  labs(title = "PCA plot: PC1 vs PC2. Including Mark's data",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2
ggplotly(fig_PC1vsPC2)

ggsave(fig_PC1vsPC2,
       file = "PCA_PC1vsPC2_PlusMark_v1.pdf",
       path = "PCA_Figures",
       width = 9, height = 6, units = "in")





