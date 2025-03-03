# Correlations between captured and not THP1 spiked samples
# E. Lamont
# 2/16/25

source("Import_data.R") # to get ProbeTest5_tpm

# Process the data a little more (not done in import data)
ProbeTest5_tpm_2 <- ProbeTest5_tpm
# Adjust the names so they are slightly shorter
names(ProbeTest5_tpm_2) <- gsub(x = names(ProbeTest5_tpm_2), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it
# add rownames to the tpm and metadata dataframes
rownames(ProbeTest5_tpm_2) <- ProbeTest5_tpm_2[,1] # add the rownames
ProbeTest5_tpm_2 <- ProbeTest5_tpm_2[,-1] # Remove the old column of rownames
# Need to make sure there is a Gene column (gets lost)
ProbeTest5_tpm_2$Gene <- rownames(ProbeTest5_tpm_2)

# Log10 transform the data
ProbeTest5_tpm_Log10 <- ProbeTest5_tpm_2 %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values


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

# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2



###################################################################
################## SINGLE GGCORRPLOT FUNCTION #####################

CAPvsNOT_ScatterCorr_Function <- function(Sample1, Sample2, my_df) {
  
  ## Sample1 goes on X-axis and is CAPTURED
  ## Sample2 goes on Y-axis and is NOT captured
  
  ScatterCorr <- my_df %>% 
    # filter(.data[[Sample2]] != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
    ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
    geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
    labs(title = paste0(Sample1, " vs ", Sample2),
         subtitle = "Pearson correlation",
         x = paste0(Sample1, " Log10(TPM+1) CAPTURED"), y = paste0(Sample2, " Log10(TPM+1) NOT captured")) + 
    stat_cor(method="pearson") + # add a correlation to the plot
    my_plot_themes
  ScatterCorr
  
}

################################################################
################## RUN FUNCTION ON SAMPLES #####################

Sample1 <- "THP1_1e6_1a" # Captured
Sample2 <- "THP1_1e6_1b" # Not Captured
my_plot <- CAPvsNOT_ScatterCorr_Function(Sample1, Sample2, ProbeTest5_tpm_Log10)
my_plot
ggsave(my_plot,
       file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")

Sample1 <- "THP1_1e6_2b" # Captured
Sample2 <- "THP1_1e6_2a" # Not Captured
my_plot <- CAPvsNOT_ScatterCorr_Function(Sample1, Sample2, ProbeTest5_tpm_Log10)
my_plot
ggsave(my_plot,
       file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")


Sample1 <- "THP1_1e6_3a" # Captured
Sample2 <- "THP1_1e6_3b" # Not Captured
my_plot <- CAPvsNOT_ScatterCorr_Function(Sample1, Sample2, ProbeTest5_tpm_Log10)
my_plot
ggsave(my_plot,
       file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")

###################################################################
####################### SINGLE GGCORRPLOT 1 #######################
Sample1 <- "THP1_1e6_1a" # Captured
Sample2 <- "THP1_1e6_1b" # Not Captured
ScatterCorr <- ProbeTest5_tpm_Log10 %>% 
  # filter(THP1_1e6_1b != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0(Sample1, " Log10(TPM+1) CAPTURED"), y = paste0(Sample2, " Log10(TPM+1) NOT captured")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")

###################################################################
####################### SINGLE GGCORRPLOT 2 #######################
Sample1 <- "THP1_1e6_2b" # Captured
Sample2 <- "THP1_1e6_2a" # Not Captured
ScatterCorr <- ProbeTest5_tpm_Log10 %>% 
  # filter(THP1_1e6_1b != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0(Sample1, " Log10(TPM+1) CAPTURED"), y = paste0(Sample2, " Log10(TPM+1) NOT captured")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")

###################################################################
####################### AVERAGES GGCORRPLOT #######################

# Take the average of the spiked capture and the averages of the spiked not capture and see if that is any better

THP1_AVERAGE_tpm_Log10 <- ProbeTest5_tpm_Log10 %>% select(Gene, THP1_1e6_1a, THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_2b, THP1_1e6_3a, THP1_1e6_3b) %>%
  mutate(
    CAPTURED_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    NOTCaptured_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b)), na.rm = TRUE),
  )

Sample1 <- "CAPTURED_THP1Spiked1e6" # Captured
Sample2 <- "NOTCaptured_THP1Spiked1e6" # Not Captured
ScatterCorr <- THP1_AVERAGE_tpm_Log10 %>% 
  # filter(NOTCaptured_THP1Spiked1e6 != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0("Log10(TPM+1) CAPTURED samples averaged"), y = paste0("Log10(TPM+1) NOT captured samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("SCALED_THP1Spiked1e6_samplesAveraged.CAPTUREDvsNOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")


###################################################################
################ AVERAGES GGCORRPLOT NOT SCALED ###################

# See what it looks like for the not scaled TPM data

# Start with ProbeTest5_tpm_NOTscaled

# Log10 transform the data
ProbeTest5_tpm_NOTscaled_Log10 <- ProbeTest5_tpm_NOTscaled %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Need to make sure there is a Gene column (gets lost)
ProbeTest5_tpm_NOTscaled_Log10$Gene <- rownames(ProbeTest5_tpm_NOTscaled_Log10)

AVERAGE_tpm_NOTscaled_Log10 <- ProbeTest5_tpm_NOTscaled_Log10 %>% select(Gene, THP1_1e6_1a, THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_2b, THP1_1e6_3a, THP1_1e6_3b) %>%
  mutate(
    CAPTURED_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    NOTCaptured_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b)), na.rm = TRUE),
  )

Sample1 <- "CAPTURED_THP1Spiked1e6" # Captured
Sample2 <- "NOTCaptured_THP1Spiked1e6" # Not Captured
ScatterCorr <- AVERAGE_tpm_NOTscaled_Log10 %>% 
  # filter(THP1_1e6_1b != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("NOT scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0("Log10(TPM+1) CAPTURED samples averaged"), y = paste0("Log10(TPM+1) NOT captured samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("NOTscaled_THP1Spiked1e6_samplesAveraged.CAPTUREDvsNOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")
