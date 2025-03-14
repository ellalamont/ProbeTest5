# Read depth for all the Sputum from ProbeTests 3,4,5
# E. Lamont 
# 2/16/25

source("Import_data.R") # to get UniqueSputum_pipeSummary

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


# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

###########################################################
############### N_GENOMIC UNIQUE SPUTUM ###################
WeekvsReads_UniqueSputum <- UniqueSputum_pipeSummary %>%
  ggplot(aes(x = Week, y = N_Genomic)) + 
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_point(aes(fill = Week, shape = Week), size = 3, alpha = 0.8, stroke = 0.5, color = "black") + # For thumbnail
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,6000000), breaks = seq(0, 6000000, 1000000)) +
  labs(title = "Unique Sputum: Week vs number reads aligned to Mtb (normal depletion)",
       subtitle = NULL, 
       x = "Weeks after start of antibiotics", 
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes
WeekvsReads_UniqueSputum

# ggplotly(WeekvsReads_sputum2)
ggsave(WeekvsReads_UniqueSputum,
       file = "WeekvsReads_UniqueSputum.pdf",
       path = "UniqueSputum_Figures",
       width = 6, height = 4, units = "in")

###########################################################
############### P_GENOMIC UNIQUE SPUTUM ###################

WeekvsPercent_UniqueSputum <- UniqueSputum_pipeSummary %>%
  ggplot(aes(x = Week, y = P_Genomic)) + 
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_point(aes(fill = Week, shape = Week), size = 3, alpha = 0.8, stroke = 0.5, color = "black") + # For thumbnail
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,75), breaks = seq(0, 75, 10)) +
  labs(title = "Unique Sputum: Week vs precent reads aligned to Mtb",
       subtitle = NULL, 
       x = "Weeks after start of antibiotics", 
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes
WeekvsPercent_UniqueSputum

# ggplotly(WeekvsReads_sputum2)
ggsave(WeekvsPercent_UniqueSputum,
       file = paste0("WeekvsPercent_UniqueSputum", ".pdf"),
       path = "UniqueSputum_Figures",
       width = 6, height = 4, units = "in")








