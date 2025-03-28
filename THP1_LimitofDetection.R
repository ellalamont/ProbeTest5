# Limit of Detection with the THP1 spiked samples graphs
# E. Lamont
# 2/18/25

source("Import_data.R") # to get LimitofDetect_pipeSummary


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
  
  
  

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

###########################################################
############# ALL CELL NUMBER VS N_GENOMIC ################

LimitofDetect_NumReads_Fig1 <- LimitofDetect_pipeSummary %>% 
  ggplot(aes(x = Ra_cells2, y = N_Genomic)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "All THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "# reads aligning to Mtb genome") + 
  # scale_y_continuous(limits = c(0,7000000), breaks = seq(0, 7000000, 1000000)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_NumReads_Fig1
ggsave(LimitofDetect_NumReads_Fig1,
       file = "All_LimitofDetect_NumReads_v1.pdf",
       path = "LimitofDetection_Figures",
       width = 6, height = 4, units = "in")

# Make it a boxplot
LimitofDetect_NumReads_Fig2 <- LimitofDetect_pipeSummary %>% 
  ggplot(aes(x = Ra_cells2, y = N_Genomic)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Run), shape = 21, alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  # geom_point(shape = 16, alpha = 0.8, size = 1.5, position = position_jitter(0.2)) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "All THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "# reads aligning to Mtb genome") + 
  scale_y_continuous(limits = c(0,17000000), breaks = seq(0, 17000000, 2000000)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_NumReads_Fig2
ggsave(LimitofDetect_NumReads_Fig2,
       file = "All_LimitofDetect_NumReads_v3.pdf",
       path = "LimitofDetection_Figures",
       width = 8, height = 4, units = "in")

###########################################################
############# ALL CELL NUMBER VS P_GENOMIC ################

LimitofDetect_PercentReads_Fig1 <- LimitofDetect_pipeSummary %>% 
  ggplot(aes(x = Ra_cells2, y = P_Genomic)) + 
  geom_point(size = 6, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  # geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "All THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "% reads aligning to Mtb genome") + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_PercentReads_Fig1
ggsave(LimitofDetect_PercentReads_Fig1,
       file = "All_LimitofDetect_PercentReads_v1.pdf",
       path = "LimitofDetection_Figures",
       width = 6, height = 4, units = "in")

# Make it a boxplot
LimitofDetect_PercentReads_Fig2 <- LimitofDetect_pipeSummary %>% 
  ggplot(aes(x = Ra_cells2, y = P_Genomic)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Run), shape = 21, alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  # geom_point(shape = 16, alpha = 0.8, size = 1.5, position = position_jitter(0.2)) + 
  labs(title = "All THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "% reads aligning to Mtb genome") + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_PercentReads_Fig2
ggsave(LimitofDetect_PercentReads_Fig2,
       file = "All_LimitofDetect_PercentReads_v3.pdf",
       path = "LimitofDetection_Figures",
       width = 8, height = 4, units = "in")


###########################################################
########## ProbeTest5 CELL NUMBER VS N_GENOMIC ############

ProbeTest5_LimitofDetect_NumReads_Fig1 <- LimitofDetect_pipeSummary %>% 
  filter(Run == "ProbeTest5") %>% 
  ggplot(aes(x = Ra_cells2, y = N_Genomic)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "ProbeTest5 THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "# reads aligning to Mtb genome") + 
  scale_y_continuous(limits = c(0,19000000), breaks = seq(0, 19000000, 2000000)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
ProbeTest5_LimitofDetect_NumReads_Fig1
ggsave(ProbeTest5_LimitofDetect_NumReads_Fig1,
       file = "ProbeTest5_LimitofDetect_NumReads_v1.pdf",
       path = "LimitofDetection_Figures",
       width = 6, height = 4, units = "in")

# ggerrorplot
ProbeTest5_LimitofDetect_NumReads_Fig2 <- LimitofDetect_pipeSummary %>% 
  filter(Run == "ProbeTest5") %>% 
  ggerrorplot(x = "Ra_cells2", y = "N_Genomic", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.8,  # Size of error bars
              add.params = list(size = 1)) +  # Size of mean points
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 5.6, y = 1000000*0.8, label = "1 million", 
           hjust = 1.1, vjust = -0.5, color = "black") + 
  # scale_y_continuous(limits = c(0,19000000), breaks = seq(0, 19000000, 2000000)) + 
  scale_y_continuous(limits = c(0,20000000), breaks = seq(0, 20000000, 4000000)) + 
  labs(# title = "ProbeTest5 THP1 cells spiked with H37Ra", 
       # subtitle = "Mean with standard deviation", 
      title = NULL, 
      subtitle = NULL, 
       x = "# bacterial cells", 
       y = "# reads aligning to \nMtb transcriptome") + 
  poster_plot_themes
ProbeTest5_LimitofDetect_NumReads_Fig2
# ggsave(ProbeTest5_LimitofDetect_NumReads_Fig2,
#        file = "ProbeTest5_LimitofDetect_NumReads_v2.pdf",
#        path = "LimitofDetection_Figures",
#        width = 7, height = 5, units = "in")
ggsave(ProbeTest5_LimitofDetect_NumReads_Fig2,
       file = "ProbeTest5_LimitofDetect_NumReads_v7.pdf",
       path = "Poster_Figures",
       width = 7.5, height = 4.5, units = "in")


###########################################################
########## ProbeTest5 CELL NUMBER VS P_GENOMIC ############

ProbeTest5_LimitofDetect_PercentReads_Fig1 <- LimitofDetect_pipeSummary %>% 
  filter(Run == "ProbeTest5") %>% 
  ggplot(aes(x = Ra_cells2, y = P_Genomic)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA) + 
  labs(title = "ProbeTest5 THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "% reads aligning to Mtb genome") + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
ProbeTest5_LimitofDetect_PercentReads_Fig1
ggsave(ProbeTest5_LimitofDetect_PercentReads_Fig1,
       file = "ProbeTest5_LimitofDetect_PercentReads_v1.pdf",
       path = "LimitofDetection_Figures",
       width = 6, height = 4, units = "in")

# ggerrorplot
ProbeTest5_LimitofDetect_PercentReads_Fig2 <- LimitofDetect_pipeSummary %>% 
  filter(Run == "ProbeTest5") %>% 
  ggerrorplot(x = "Ra_cells2", y = "P_Genomic", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black") + 
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  labs(title = "ProbeTest5 THP1 cells spiked with H37Ra", 
       subtitle = "Mean with standard deviation", 
       x = "# spiked in H37Ra cells", 
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes
ProbeTest5_LimitofDetect_PercentReads_Fig2
ggsave(ProbeTest5_LimitofDetect_PercentReads_Fig2,
       file = "ProbeTest5_LimitofDetect_PercentReads_v2.pdf",
       path = "LimitofDetection_Figures",
       width = 7, height = 5, units = "in")



###########################################################
########## ALL CELL NUMBER VS AtLeast.10.Reads ############

# Make it a boxplot
LimitofDetect_10Reads_Fig1 <- LimitofDetect_pipeSummary %>% 
  ggplot(aes(x = Ra_cells2, y = AtLeast.10.Reads)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Run), shape = 21, alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  # geom_point(shape = 16, alpha = 0.8, size = 1.5, position = position_jitter(0.2)) + 
  geom_hline(yintercept = 4499*0.8, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 0.9, y = 4499*0.8, label = "80%", 
           hjust = 1.1, vjust = -0.5, color = "black") + 
  geom_hline(yintercept = 4499*0.5, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 0.9, y = 4499*0.5, label = "50%", 
           hjust = 1.1, vjust = -0.5, color = "black") + 
  labs(title = "All THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "# genes with at least 10 reads aligning") + 
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, 500)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
LimitofDetect_10Reads_Fig1
ggsave(LimitofDetect_10Reads_Fig1,
       file = "All_LimitofDetect_AtLeast10Reads_v1.pdf",
       path = "LimitofDetection_Figures",
       width = 8, height = 4, units = "in")

###########################################################
####### ProbeTest5 CELL NUMBER VS AtLeast.10.Reads ########

ProbeTest5_LimitofDetect_10Reads_Fig1 <- LimitofDetect_pipeSummary %>% 
  filter(Run == "ProbeTest5") %>%
  ggplot(aes(x = Ra_cells2, y = AtLeast.10.Reads)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, fill = "#03A9F4", shape = 21) + 
  geom_hline(yintercept = 4499*0.8, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 0.9, y = 4499*0.8, label = "80%", 
           hjust = 1.1, vjust = -0.5, color = "black") + 
  geom_hline(yintercept = 4499*0.5, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 0.9, y = 4499*0.5, label = "50%", 
           hjust = 1.1, vjust = -0.5, color = "black") + 
  labs(title = "ProbeTest5 THP1 cells spiked with H37Ra",
       subtitle = NULL, 
       x = "# spiked in H37Ra cells", 
       y = "# genes with at least 10 reads aligning") + 
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, 500)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
ProbeTest5_LimitofDetect_10Reads_Fig1
ggsave(ProbeTest5_LimitofDetect_10Reads_Fig1,
       file = "ProbeTest5_LimitofDetect_AtLeast10Reads_v1.pdf",
       path = "LimitofDetection_Figures",
       width = 6, height = 4, units = "in")


# ggerrorplot
ProbeTest5_LimitofDetect_10Reads_Fig2 <- LimitofDetect_pipeSummary %>% 
  filter(Run == "ProbeTest5") %>% 
  ggerrorplot(x = "Ra_cells2", y = "AtLeast.10.Reads", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black") + 
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, 500)) + 
  geom_hline(yintercept = 4499*0.8, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 0.9, y = 4499*0.8, label = "80%", 
           hjust = 1.1, vjust = -0.5, color = "black") + 
  geom_hline(yintercept = 4499*0.5, linetype = "dashed", alpha = 0.5) + 
  annotate("text", x = 0.9, y = 4499*0.5, label = "50%", 
           hjust = 1.1, vjust = -0.5, color = "black") + 
  labs(title = "ProbeTest5 THP1 cells spiked with H37Ra", 
       subtitle = "Mean with standard deviation", 
       x = "# spiked in H37Ra cells", 
       y = "# genes with at least 10 reads aligning") + 
  my_plot_themes
ProbeTest5_LimitofDetect_10Reads_Fig2
ggsave(ProbeTest5_LimitofDetect_10Reads_Fig2,
       file = "ProbeTest5_LimitofDetect_AtLeast10Reads_v2.pdf",
       path = "LimitofDetection_Figures",
       width = 7, height = 5, units = "in")
