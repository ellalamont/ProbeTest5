# Looking at the THP1 spiked samples, how many genes have at least 10 or 100 reads
# E. Lamont
# 3/7/25

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

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

###########################################################
#################### ProbeTest5 TPM #######################

ProbeTest5_tpm_2 <- ProbeTest5_tpm
# Adjust the names so they are slightly shorter
names(ProbeTest5_tpm_2) <- gsub(x = names(ProbeTest5_tpm_2), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

CapturedVsNot_SampleIDs <- c("THP1_1e6_1a", "THP1_1e6_1b", "THP1_1e6_2a", "THP1_1e6_2b", "THP1_1e6_3a"  , "THP1_1e6_3b")

CapturedVsNot_tpm <- ProbeTest5_tpm_2 %>% select(all_of(CapturedVsNot_SampleIDs))

###########################################################
################# HOW MANY HAVE >1 READ ###################

# See how many have at least one read

counts <- CapturedVsNot_tpm %>%
  summarise(across(everything(), ~ sum(. > 0)))
result <- 4499-counts # This is how many genes have 0 reads aligning in each!
# THP1_1e6_1a THP1_1e6_1b THP1_1e6_2a THP1_1e6_2b THP1_1e6_3a THP1_1e6_3b
# 1          14         265         247          15          12         923


###########################################################
########## ALL CELL NUMBER VS AtLeast.10.Reads ############
# 4/22/25

#ggerrroplot
ProbeTest5_CapturedvsNot_10Reads_Fig1 <- CapturedVsNot_pipeSummary %>% 
  # filter(Run == "ProbeTest5") %>% 
  ggerrorplot(x = "Probe", y = "AtLeast.10.Reads", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black") + 
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
  scale_x_discrete(labels = c("None" = "Not captured",
                              "JA2" = "Captured")) + 
  my_plot_themes
ProbeTest5_CapturedvsNot_10Reads_Fig1
ggsave(ProbeTest5_CapturedvsNot_10Reads_Fig1,
       file = "AtLeast10Reads_CapturedVsNot.pdf",
       path = "CaptureVsNot_Figures",
       width = 6, height = 4, units = "in")



