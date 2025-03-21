# Correlations between captured and not THP1 spiked samples
# E. Lamont
# 2/16/25

source("Import_data.R") # to get THP1Spike_tpm_CorrectScales

# Want the captured samples to be scaled and the not captured samples to be not scaled
# Captured = "THP1_1e6_1a", "THP1_1e6_2b", "THP1_1e6_3a" 
# Not Captured = "THP1_1e6_1b", "THP1_1e6_2a", "THP1_1e6_3b" 

# Log10 transform the data
THP1Spike_tpm_CorrectScales_Log10 <- THP1Spike_tpm_CorrectScales %>% 
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




# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2



###################################################################
################## SINGLE GGCORRPLOT FUNCTION #####################

CAPvsNOT_ScatterCorr_Function <- function(Sample1, Sample2, my_df) {
  
  ## Sample1 goes on X-axis and is CAPTURED & Scaled
  ## Sample2 goes on Y-axis and is NOT captured & NOT scaled
  
  ScatterCorr <- my_df %>% 
    # filter(.data[[Sample2]] != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
    ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
    geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
    labs(title = paste0(Sample1, " vs ", Sample2),
         subtitle = "Pearson correlation",
         x = paste0(Sample1, " Log10(TPM+1) CAPTURED & scaled"), y = paste0(Sample2, " Log10(TPM+1) NOT captured & NOT scaled")) + 
    stat_cor(method="pearson") + # add a correlation to the plot
    my_plot_themes
  ScatterCorr
  
}

################################################################
################## RUN FUNCTION ON SAMPLES #####################

Sample1 <- "THP1_1e6_1a" # Captured
Sample2 <- "THP1_1e6_1b" # Not Captured
my_plot <- CAPvsNOT_ScatterCorr_Function(Sample1, Sample2, THP1Spike_tpm_CorrectScales_Log10)
my_plot
ggsave(my_plot,
       file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations_CorrectScales",
       width = 7, height = 5, units = "in")

Sample1 <- "THP1_1e6_2b" # Captured
Sample2 <- "THP1_1e6_2a" # Not Captured
my_plot <- CAPvsNOT_ScatterCorr_Function(Sample1, Sample2, THP1Spike_tpm_CorrectScales_Log10)
my_plot
ggsave(my_plot,
       file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations_CorrectScales",
       width = 7, height = 5, units = "in")


Sample1 <- "THP1_1e6_3a" # Captured
Sample2 <- "THP1_1e6_3b" # Not Captured
my_plot <- CAPvsNOT_ScatterCorr_Function(Sample1, Sample2, THP1Spike_tpm_CorrectScales_Log10)
my_plot
ggsave(my_plot,
       file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations_CorrectScales",
       width = 7, height = 5, units = "in")

###################################################################
####################### SINGLE GGCORRPLOT 1 #######################
# Sample1 <- "THP1_1e6_1a" # Captured
# Sample2 <- "THP1_1e6_1b" # Not Captured
# ScatterCorr <- ProbeTest5_tpm_Log10 %>% 
#   # filter(THP1_1e6_1b != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
#   ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
#   geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
#   geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
#   labs(title = paste0(Sample1, " vs ", Sample2),
#        subtitle = "Pearson correlation",
#        x = paste0(Sample1, " Log10(TPM+1) CAPTURED"), y = paste0(Sample2, " Log10(TPM+1) NOT captured")) + 
#   stat_cor(method="pearson") + # add a correlation to the plot
#   my_plot_themes
# ScatterCorr
# # ggplotly(ScatterCorr)
# ggsave(ScatterCorr,
#        file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
#        path = "Correlation_Figures/THP1Spiked_Correlations",
#        width = 7, height = 5, units = "in")
# 
# ###################################################################
# ####################### SINGLE GGCORRPLOT 2 #######################
# Sample1 <- "THP1_1e6_2b" # Captured
# Sample2 <- "THP1_1e6_2a" # Not Captured
# ScatterCorr <- ProbeTest5_tpm_Log10 %>% 
#   # filter(THP1_1e6_1b != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
#   ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
#   geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
#   geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
#   labs(title = paste0(Sample1, " vs ", Sample2),
#        subtitle = "Pearson correlation",
#        x = paste0(Sample1, " Log10(TPM+1) CAPTURED"), y = paste0(Sample2, " Log10(TPM+1) NOT captured")) + 
#   stat_cor(method="pearson") + # add a correlation to the plot
#   my_plot_themes
# ScatterCorr
# # ggplotly(ScatterCorr)
# ggsave(ScatterCorr,
#        file = paste0(Sample1, ".CAPTURED_ComparedTo_", Sample2, ".NOTCaptured.pdf"),
#        path = "Correlation_Figures/THP1Spiked_Correlations",
#        width = 7, height = 5, units = "in")

###################################################################
####################### AVERAGES GGCORRPLOT #######################

# Take the average of the spiked capture and the averages of the spiked not capture and see if that is any better

THP1Spike_AVERAGE_tpm_CorrectScales_Log10 <- THP1Spike_tpm_CorrectScales_Log10 %>% select(Gene, THP1_1e6_1a, THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_2b, THP1_1e6_3a, THP1_1e6_3b) %>%
  mutate(
    CAPTURED_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    NOTCaptured_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b)), na.rm = TRUE),
  )

Sample1 <- "CAPTURED_THP1Spiked1e6" # Captured
Sample2 <- "NOTCaptured_THP1Spiked1e6" # Not Captured
ScatterCorr <- THP1Spike_AVERAGE_tpm_CorrectScales_Log10 %>% 
  # filter(NOTCaptured_THP1Spiked1e6 != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; Captured Scaled; Not captured Not scaled",
       x = paste0("Log10(TPM+1) CAPTURED samples averaged"), y = paste0("Log10(TPM+1) NOT captured samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("SCALED_THP1Spiked1e6_samplesAveraged.CAPTUREDvsNOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations_CorrectScales",
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

# For Poster
ScatterCorr_poster <- AVERAGE_tpm_NOTscaled_Log10 %>% 
  # filter(THP1_1e6_1b != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.55, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = NULL,
       subtitle = NULL,
       x = paste0("Captured samples Log10(TPM+1)"), y = paste0("Not captured samples \nLog10(TPM+1)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  poster_plot_themes
ScatterCorr_poster
# ggplotly(ScatterCorr)
ggsave(ScatterCorr_poster,
       file = paste0("NOTscaled_THP1Spiked1e6_samplesAveraged.CAPTUREDvsNOTCaptured_v2.pdf"),
       path = "Poster_Figures",
       width = 7.5, height = 4.5, units = "in")






# Get the individual plots
Sample1 <- "THP1_1e6_1a" # Captured
Sample2 <- "THP1_1e6_1b" # Not Captured
Sample1 <- "THP1_1e6_2b" # Captured
Sample2 <- "THP1_1e6_2a" # Not Captured
Sample1 <- "THP1_1e6_3a" # Captured
Sample2 <- "THP1_1e6_3b" # Not Captured
ScatterCorr <- AVERAGE_tpm_NOTscaled_Log10 %>% 
  # filter(THP1_1e6_1b != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("NOT scaled Samples: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0("Log10(TPM+1) CAPTURED ", Sample1), y = paste0("Log10(TPM+1) NOT captured ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
# ggsave(ScatterCorr,
#        file = paste0("NOTscaled_", Sample1, "vs", Sample2, ".pdf"),
#        path = "Correlation_Figures/THP1Spiked_Correlations",
#        width = 7, height = 5, units = "in")






###################################################################
############# AVERAGES GGCORRPLOT NOT logTransformed ##############

# Scaled data but not log transformed

# Start with ProbeTest5_tpm_2

AVERAGE_tpm_scaled <- ProbeTest5_tpm_2 %>% select(Gene, THP1_1e6_1a, THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_2b, THP1_1e6_3a, THP1_1e6_3b) %>%
  mutate(
    CAPTURED_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    NOTCaptured_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b)), na.rm = TRUE),
  )

Sample1 <- "CAPTURED_THP1Spiked1e6" # Captured
Sample2 <- "NOTCaptured_THP1Spiked1e6" # Not Captured
ScatterCorr <- AVERAGE_tpm_scaled %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("NOT Transformed, Yes scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0("TPM CAPTURED samples averaged"), y = paste0("TPM NOT captured samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  scale_x_continuous(limits = c(0,20000), breaks = seq(0, 20000, 4000)) + 
  scale_y_continuous(limits = c(0,20000), breaks = seq(0, 20000, 4000)) + 
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("NOT.logTrans_THP1Spiked1e6_samplesAveraged.CAPTUREDvsNOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")


#### See what the individual comparisons look like...

Sample1 <- "THP1_1e6_1a" # Captured
Sample2 <- "THP1_1e6_1b" # Not Captured
ScatterCorr <- ProbeTest5_tpm_2 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("NOT Transformed, Yes scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0(Sample1, " TPM CAPTURED"), y = paste0(Sample2, " TPM NOT captured")) +
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr

Sample1 <- "THP1_1e6_2b" # Captured
Sample2 <- "THP1_1e6_2a" # Not Captured
ScatterCorr <- ProbeTest5_tpm_2 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("NOT Transformed, Yes scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0(Sample1, " TPM CAPTURED"), y = paste0(Sample2, " TPM NOT captured")) +
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr


Sample1 <- "THP1_1e6_3a" # Captured
Sample2 <- "THP1_1e6_3b" # Not Captured
ScatterCorr <- ProbeTest5_tpm_2 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("NOT Transformed, Yes scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0(Sample1, " TPM CAPTURED"), y = paste0(Sample2, " TPM NOT captured")) +
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr



###################################################################
####### AVERAGES GGCORRPLOT NOT logTransformed NOT Scaled #########

# Start with ProbeTest5_tpm_NOTscaled
ProbeTest5_tpm_NOTscaled$Gene <- rownames(ProbeTest5_tpm_NOTscaled)

AVERAGE_tpm_NOTscaled <- ProbeTest5_tpm_NOTscaled %>% select(Gene, THP1_1e6_1a, THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_2b, THP1_1e6_3a, THP1_1e6_3b) %>%
  mutate(
    CAPTURED_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    NOTCaptured_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b)), na.rm = TRUE),
  )

Sample1 <- "CAPTURED_THP1Spiked1e6" # Captured
Sample2 <- "NOTCaptured_THP1Spiked1e6" # Not Captured
ScatterCorr <- AVERAGE_tpm_NOTscaled %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("NOT Transformed, NOT scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0("TPM CAPTURED samples averaged"), y = paste0("TPM NOT captured samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("NOTscaled_NOT.logTrans_THP1Spiked1e6_samplesAveraged.CAPTUREDvsNOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")






###################################################################
############### AVERAGES GGCORRPLOT LOG2 TRANSFORMED ##############

# Scaled data 
# Betin does a log2 transformation so I just want to see how mine compares

# Start with ProbeTest5_tpm_2
# Log2 transform the data
ProbeTest5_tpm_Log2 <- ProbeTest5_tpm_2 %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log2(.x))) # Log transform the values


AVERAGE_tpm_Log2_scaled <- ProbeTest5_tpm_Log2 %>% select(Gene, THP1_1e6_1a, THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_2b, THP1_1e6_3a, THP1_1e6_3b) %>%
  mutate(
    CAPTURED_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    NOTCaptured_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b)), na.rm = TRUE),
  )

Sample1 <- "CAPTURED_THP1Spiked1e6" # Captured
Sample2 <- "NOTCaptured_THP1Spiked1e6" # Not Captured
ScatterCorr <- AVERAGE_tpm_Log2_scaled %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("NOT Transformed, Yes scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0("Log2(TPM+1) CAPTURED samples averaged"), y = paste0("Log2(TPM+1) NOT captured samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,20000), breaks = seq(0, 20000, 4000)) + 
  # scale_y_continuous(limits = c(0,20000), breaks = seq(0, 20000, 4000)) + 
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("Log2_THP1Spiked1e6_samplesAveraged.CAPTUREDvsNOTCaptured.pdf"),
       path = "Correlation_Figures/THP1Spiked_Correlations",
       width = 7, height = 5, units = "in")







