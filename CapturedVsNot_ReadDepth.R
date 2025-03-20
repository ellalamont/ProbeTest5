# Read depth for all captured vs not THP1 samples for ProbeTest5
# E. Lamont 
# 2/24/25

source("Import_data.R") # to get CapturedVsNot_pipeSummary

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none", legend.text=element_text(size=10),
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

poster_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=20), 
        axis.text.x = element_text(angle = 0, size=20, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=18),
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
###################### N_GENOMIC ##########################
WeekvsReads_CapturedVsNot <- CapturedVsNot_pipeSummary %>%
  ggplot(aes(x = Probe, y = N_Genomic)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Replicates), color = "black", size = 0.5, linetype = "dashed") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,19000000), breaks = seq(0, 19000000, 2000000)) +
  labs(title = "THP1 spiked with 1e6 H37Ra",
       subtitle = NULL, 
       x = NULL, 
       y = "# reads aligning to Mtb genome") + 
  scale_x_discrete(labels = c("None" = "Not captured",
                              "JA2" = "Captured")) + 
  my_plot_themes
WeekvsReads_CapturedVsNot

# ggplotly(WeekvsReads_sputum2)
ggsave(WeekvsReads_CapturedVsNot,
       file = "WeekvsReads_CapturedVsNot.pdf",
       path = "CaptureVsNot_Figures",
       width = 6, height = 4, units = "in")


# For poster
WeekvsReads_CapturedVsNot_poster <- CapturedVsNot_pipeSummary %>%
  ggplot(aes(x = Probe, y = N_Genomic)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Replicates), color = "black", size = 0.5, linetype = "dashed") + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,20000000), breaks = seq(0, 20000000, 4000000)) +
  labs(title = "THP1 spiked with 1e6 H37Ra",
       subtitle = "subtitle to remove", 
       x = "x axis label to remove", 
       y = "# reads aligning to Mtb transcriptome") + 
  scale_x_discrete(labels = c("None" = "Not captured",
                              "JA2" = "Transcript captured")) + 
  poster_plot_themes
WeekvsReads_CapturedVsNot_poster
# ggplotly(WeekvsReads_sputum2)
ggsave(WeekvsReads_CapturedVsNot_poster,
       file = "WeekvsReads_CapturedVsNot_poster.pdf",
       path = "Poster_Figures",
       width = 7, height = 4.5, units = "in")


###########################################################
###################### P_GENOMIC ##########################
WeekvsPercent_CapturedVsNot <- CapturedVsNot_pipeSummary %>%
  ggplot(aes(x = Probe, y = P_Genomic)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Replicates), color = "black", size = 0.5, linetype = "dashed") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  labs(title = "THP1 spiked with 1e6 H37Ra",
       subtitle = NULL, 
       x = NULL, 
       y = "% reads aligning to Mtb genome") + 
  scale_x_discrete(labels = c("None" = "Not captured",
                              "JA2" = "Captured")) + 
  my_plot_themes
WeekvsPercent_CapturedVsNot

# ggplotly(WeekvsReads_sputum2)
ggsave(WeekvsPercent_CapturedVsNot,
       file = "WeekvsPercent_CapturedVsNot.pdf",
       path = "CaptureVsNot_Figures",
       width = 6, height = 4, units = "in")






