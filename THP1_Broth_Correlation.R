# Compare captured THP1 to uncaptured broth (all not scaled)
# E. Lamont 
# 3/14/25

source("Import_data.R") # to get ProbeTest5_tpm_NOTscaled


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
        # plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )

poster_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
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


my_sampleIDs <- c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6",
                  "THP1_1e6_1a", "THP1_1e6_2b", "THP1_1e6_3a", "Gene")

my_tpm_subset <- ProbeTest5_tpm_NOTscaled %>% select(all_of(my_sampleIDs))

# Log10 transform the data
my_tpm_subset_Log10 <- my_tpm_subset %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values


# Add columns for averages
my_tpm_subset_Log10 <- my_tpm_subset_Log10 %>% 
  mutate(
    AVERAGE_THP1Spiked = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    AVERAGE_BrothNotCaptured = rowMeans(select(., c(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6)), na.rm = TRUE),
  )


###################################################################
####################### AVERAGES GGCORRPLOT #######################

Sample1 <- "AVERAGE_THP1Spiked" # THP1 spiked Captured
Sample2 <- "AVERAGE_BrothNotCaptured" # Broth Not Captured
ScatterCorr <- my_tpm_subset_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; 1e6 Ra THP1 spiked captured VS Broth Not captured (Not scaled)",
       x = paste0("Log10(TPM+1) Ra1e6 THP1 samples averaged"), y = paste0("Log10(TPM+1) NOT captured Broth samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("THP1Spiked1e6_vs_BrothNOTCaptured_Averages.pdf"),
       path = "Correlation_Figures",
       width = 7, height = 5, units = "in")

# For poster
ScatterCorr_poster <- my_tpm_subset_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.55, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = NULL,
       subtitle = NULL,
       x = paste0("Captured samples Log10(TPM+1)"), y = paste0("Not captured samples (broth) \nLog10(TPM+1)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  poster_plot_themes
ScatterCorr_poster
# ggplotly(ScatterCorr)
ggsave(ScatterCorr_poster,
       file = paste0("THP1Spiked1e6_vs_BrothNOTCaptured_Averages.pdf"),
       path = "Poster_Figures",
       width = 7.5, height = 4.5, units = "in")



