# Make Pie chart of where reads are going 
# E. Lamont
# 2/24/25

# https://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization
# https://r-graph-gallery.com/pie-plot.html
# https://ggplot2.tidyverse.org/reference/coord_polar.html



source("Import_data.R") # to get CapturedVsNot_pipeSummary

# Plot basics
my_plot_themes <- theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        strip.text = element_text(size = 12, face = "bold"), # For the facet
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10, hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_blank(), # Remove facet panel borders
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank())


###########################################################
###################### REORDER DF #########################

my_pipeSummary_ReadCount <- CapturedVsNot_pipeSummary %>% 
  mutate(N_NotMtb = (RawReads - N_RiboClear - N_Genomic)) %>% 
  pivot_longer(cols = c("N_RiboClear", "N_Genomic", "N_NotMtb"), names_to = "Number", values_to = "ReadCount") %>% 
  mutate(ReadPercent = round((ReadCount/RawReads)*100, 1)) 
# Why didn't I just use P??

replacement_values <- c(N_Genomic = "mRNA", N_RiboClear = "rRNA", N_NotMtb = "other RNA")
my_pipeSummary_ReadCount <- my_pipeSummary_ReadCount %>% 
  mutate(Number = replacement_values[Number])




###########################################################
##################### MULTIPLE PIE ########################

PieChart_v1 <- my_pipeSummary_ReadCount %>% 
  filter(SampleID %in% c("THP1_1e6_1a", "THP1_1e6_1b")) %>%
  mutate(SampleID = factor(SampleID, levels = c("THP1_1e6_1b", "THP1_1e6_1a"))) %>% 
  arrange(SampleID, desc(Number)) %>%  # Need this so the numbers go to the correct slices
  group_by(SampleID) %>% # Need this or it won't be a round pie! 
  mutate(cumulative = cumsum(ReadPercent), midpoint = cumulative - ReadPercent / 2) %>%
  ungroup() %>%
  ggplot(aes(x = "", y = ReadPercent, fill = Number)) +
  geom_bar(width = 1, stat = "identity", color = "black") + 
  coord_polar(theta = "y", start = 0) + 
  facet_wrap(~SampleID, labeller = as_labeller(c("THP1_1e6_1b" = "Not captured", "THP1_1e6_1a" = "Captured"))) +
  scale_fill_manual(values = c("#00CED1", "#708090", "#E0D8B0")) + 
  # geom_text(aes(y = midpoint, label = paste(Number, "\n", scales::percent(ReadPercent / 100))), size = 4, color = "black") +
  geom_text_repel(aes(y = midpoint, label = paste(Number, "\n", scales::percent(ReadPercent / 100))), size = 4, color = "black", box.padding = 0.3, force = 2, force_pull = 2, min.segment.length = 0.2, segment.size = 0.5) + 
  labs(title = "THP1 cells spiked with 1e6 H37Ra") + 
  my_plot_themes
PieChart_v1
ggsave(PieChart_v1,
       file = "THP1_1e6_1_CaptureVsNone_PIE.pdf",
       path = "Pie_Figures",
       width = 6, height = 4, units = "in")




###########################################################
################## AVERAGE REPLICATES #####################

Averages_pipeSummary <- CapturedVsNot_pipeSummary %>%
  group_by(Probe) %>%
  summarize(count = n(),
            mean_P_Genomic = mean(P_Genomic),
            mean_P_RiboClear = mean(P_RiboClear),
            mean_P_NoHit = mean(P_NoHit)) %>%
  pivot_longer(cols = c("mean_P_Genomic", "mean_P_RiboClear", "mean_P_NoHit"), names_to = "Percent_Type", values_to = "Percent") %>%
  mutate(Percent = round(Percent, 1)) 

replacement_values <- c(mean_P_Genomic = "mRNA", mean_P_RiboClear = "rRNA", mean_P_NoHit = "other RNA")
Averages_pipeSummary <- Averages_pipeSummary %>% 
  mutate(Percent_Type = replacement_values[Percent_Type])



###########################################################
############## MULTIPLE PIE - AVERAGES ####################

PieChart_Averages <- Averages_pipeSummary %>% 
  filter(Probe %in% c("JA2", "None")) %>%
  mutate(Probe = factor(Probe, levels = c("None", "JA2"))) %>% 
  arrange(Probe, desc(Percent_Type)) %>%  # Need this so the numbers go to the correct slices
  group_by(Probe) %>% # Need this or it won't be a round pie! 
  mutate(cumulative = cumsum(Percent), midpoint = cumulative - Percent / 2) %>%
  ungroup() %>%
  ggplot(aes(x = "", y = Percent, fill = Percent_Type)) +
  geom_bar(width = 1, stat = "identity", color = "black") + 
  coord_polar(theta = "y", start = 0) + 
  facet_wrap(~Probe, labeller = as_labeller(c("None" = "Not captured average", "JA2" = "Captured average"))) +
  scale_fill_manual(values = c("#00CED1", "#708090", "#E0D8B0")) + 
  geom_text_repel(aes(y = midpoint, label = paste(Percent_Type, "\n", scales::percent(Percent / 100))), size = 4, color = "black", box.padding = 0.3, force = 2, force_pull = 2, min.segment.length = 0.2, segment.size = 0.5) + 
  labs(title = "AVERAGES THP1 cells spiked with 1e6 H37Ra") + 
  my_plot_themes
PieChart_Averages
ggsave(PieChart_Averages,
       file = "Averages_CaptureVsNone_PIE.pdf",
       path = "Pie_Figures",
       width = 6, height = 4, units = "in")

