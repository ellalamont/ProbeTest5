# Correlations between captured and not captured Broth samples
# E. Lamont
# 2/16/25

source("Import_data.R") # to get Broth_tpm
Broth_tpm$Gene <- rownames(Broth_tpm)
Broth_tpm_NOTscaled$Gene <- rownames(Broth_tpm_NOTscaled)

# Log10 transform the data
Broth_tpm_Log10 <- Broth_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values
Broth_tpm_NOTscaled_Log10 <- Broth_tpm_NOTscaled %>% 
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

###########################################################
############ LOOP THROUGH ALL SCALED SAMPLES ##############

my_path <- "Correlation_Figures/Broth_Correlations"

for (i in 1:(length(colnames(Broth_tpm_Log10)) -1 -1)) {
  for (j in (i + 1):(length(colnames(Broth_tpm_Log10))-1)) {
    if (i != j) { # Avoid comparing the same samples or repeating comparisons
    # Access the samples
      Sample1 <- colnames(Broth_tpm_Log10[i])
      Sample2 <- colnames(Broth_tpm_Log10[j])
      # cat("Comparing:", Sample1, "with", Sample2, "\n")
      filename <- paste0(Sample1, "_ComparedTo_", Sample2, ".pdf")

      ScatterCorr <- Broth_tpm_Log10 %>%
        ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) +
        geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
        geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
        labs(title = paste0(Sample1, " vs ", Sample2),
             subtitle = "Pearson correlation",
             x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) +
        stat_cor(method="pearson") + # add a correlation to the plot
        my_plot_themes

      ggsave(ScatterCorr,
             file = filename,
             path = my_path,
             width = 7, height = 5, units = "in")
    }
  }
}
# This works but it obviously making a bunch that I don't need, will go through and delete, so don't rerun!!

###########################################################
########### LOOP THROUGH ALL NOT scaled SAMPLES ###########

my_path <- "Correlation_Figures/Broth_Correlations_NOTscaled"

for (i in 1:(length(colnames(Broth_tpm_NOTscaled_Log10)) -1 -1)) {
  for (j in (i + 1):(length(colnames(Broth_tpm_NOTscaled_Log10))-1)) {
    if (i != j) { # Avoid comparing the same samples or repeating comparisons
      # Access the samples
      Sample1 <- colnames(Broth_tpm_NOTscaled_Log10[i])
      Sample2 <- colnames(Broth_tpm_NOTscaled_Log10[j])
      # cat("Comparing:", Sample1, "with", Sample2, "\n")
      filename <- paste0(Sample1, "_ComparedTo_", Sample2, "_NOTscaled.pdf")
      
      ScatterCorr <- Broth_tpm_NOTscaled_Log10 %>%
        ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) +
        geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
        geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
        labs(title = paste0(Sample1, " vs ", Sample2, " NOT scaled TPM"),
             subtitle = "Pearson correlation",
             x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) +
        stat_cor(method="pearson") + # add a correlation to the plot
        my_plot_themes
      
      ggsave(ScatterCorr,
             file = filename,
             path = my_path,
             width = 7, height = 5, units = "in")
    }
  }
}


###################################################################
######################### SINGLE GGCORRPLOT #######################
Sample1 <- "H37Ra_Broth_3" # Captured
Sample2 <- "H37Ra_Broth_4" # Not Captured
ScatterCorr <- Broth_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0(Sample1, " Log10(TPM) CAPTURED"), y = paste0(Sample2, " Log10(TPM) NOT captured")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0(Sample1, "_ComparedTo_", Sample2, "_v2.pdf"),
       path = "Correlation_Figures/Broth_Correlations",
       width = 7, height = 5, units = "in")


ScatterCorr <- Broth_tpm %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0(Sample1, " CAPTURED"), y = paste0(Sample2, " NOT captured")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0(Sample1, "_ComparedTo_", Sample2, "_notTransformed.pdf"),
       path = my_path,
       width = 7, height = 5, units = "in")



###################################################################
############ Loess (Local Regression) Normalization ###############
Sample1 <- "H37Ra_Broth_3" # Captured
Sample2 <- "H37Ra_Broth_4" # Not Captured

fit <- loess(H37Ra_Broth_3 ~ H37Ra_Broth_4, data = Broth_tpm_Log10)
# Predict the trend
predicted <- predict(fit, newdata = Broth_tpm_Log10$H37Ra_Broth_3)
# Calculate residuals (or adjusted values)
Broth_tpm_Log10$H37Ra_Broth_3_adjusted <- Broth_tpm_Log10$H37Ra_Broth_3 - (predicted - mean(Broth_tpm_Log10$H37Ra_Broth_3))

ScatterCorr <- Broth_tpm_Log10 %>% 
  ggplot(aes(x = H37Ra_Broth_3_adjusted, y = H37Ra_Broth_4)) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0(Sample1, " Log10(TPM) CAPTURED"), y = paste0(Sample2, " Log10(TPM) NOT captured")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr



