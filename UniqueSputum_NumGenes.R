# Looking at the sputum samples, how many genes have at least 10 or 100 reads
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
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default


###########################################################
########### AT LEAST 10 GENES UNIQUE SPUTUM ###############

UniqueSputum_10Genes <- UniqueSputum_pipeSummary %>% 
  ggplot(aes(x = Week, y = AtLeast.10.Reads, text = SampleID)) + 
  geom_point(aes(fill = Week, shape = Week), size = 6, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_text_repel(aes(label = format(AtLeast.10.Reads, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_fill_manual(values=c(`0` = "#0072B2", `2` = "#E66900", `4`= "#009E73")) +  
  scale_shape_manual(values=c(`0` = 21, `2` = 22, `4`= 23)) + 
  scale_y_continuous(limits = c(0,4499), breaks = seq(0, 4500, 500)) + 
  labs(title = "Unique Sputum # genes aligned to Mtb with at least 10 TPM (dotted line=80%)",
       subtitle = NULL,
       x = "Weeks after start of antibiotics", 
       y = "# genes with at least 10 reads") + 
  geom_hline(yintercept = 4499*0.8, linetype = "dashed", alpha = 0.5) + 
  my_plot_themes
UniqueSputum_10Genes

ggsave(UniqueSputum_10Genes,
       file = "UniqueSputum_10Genes.pdf",
       path = "UniqueSputum_Figures",
       width = 6, height = 4, units = "in")
ggsave(UniqueSputum_10Genes,
       file = "UniqueSputum_10Genes.png",
       path = "UniqueSputum_Figures",
       width = 6, height = 4, units = "in")





