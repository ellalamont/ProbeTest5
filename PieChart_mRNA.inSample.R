# Making a nice pie chart to visualize how much bacterial mRNA is in a sample
# E. Lamont 
# 3/5/25

# Using data from 
# https://pubs.acs.org/doi/10.1021/acsinfecdis.8b00369
# Basically re-creating the pie chart in this paper
# Another paper is: https://www.nature.com/articles/nrmicro2852?utm_source=acs&getft_integrator=acs


# Plot basics
my_plot_themes <- theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"), # For the facet
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 10, hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        plot.subtitle = element_text(size = 10), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.border = element_blank(),  # Remove facet panel borders
        legend.background = element_blank(),  # Ensure the legend background is blank
        legend.box.background = element_blank())  # Ensure the legend box has no border


###########################################################
#################### MAKE THE DATA ########################

# Using data from 
# https://pubs.acs.org/doi/10.1021/acsinfecdis.8b00369

euk_data <- data.frame(
  RNA_Type = factor(c("eukaryotic rRNA", "eukaryotic tRNA/ncRNA/etc", "eukaryotic mRNA", "bacterial total RNA")),
  Percent = c(80, 12, 7, 1)
)

bact_data <- data.frame(
  RNA_Type = c("bacterial rRNA", "bacterial tRNA/ncRNA/etc", "bacterial mRNA"),
  Percent = c(85, 10, 5)
)

###########################################################
#################### EUKARYOTIC PIE #######################


# Define custom colors for RNA_Type
rna_colors <- c("eukaryotic rRNA" = "#707070", 
                "eukaryotic tRNA/ncRNA/etc" = "#A0A0A0", 
                "eukaryotic mRNA" = "#E0E0E0", 
                "bacterial total RNA" = "#D95319")

euk_data <- data.frame(
  RNA_Type = factor(c("eukaryotic rRNA", "bacterial total RNA", "eukaryotic tRNA/ncRNA/etc", 
                      "eukaryotic mRNA"), 
                    levels = c("eukaryotic rRNA", "bacterial total RNA", "eukaryotic tRNA/ncRNA/etc", "eukaryotic mRNA")),
  Percent = c(80, 1, 12, 7)
)

PieChart_euk <- euk_data %>% 
  arrange(RNA_Type, desc(Percent)) %>%  # Need this so the numbers go to the correct slices
  group_by(RNA_Type) %>% # Need this or it won't be a round pie! 
  ungroup() %>%
  ggplot(aes(x = "", y = Percent, fill = RNA_Type)) +
  geom_bar(width = 1, stat = "identity", color = "black") + 
  coord_polar(theta = "y", start = 0) + 
  scale_fill_manual(values = rna_colors) + 
  geom_text(aes(label = paste0(Percent, "%")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black") + 
  labs(title = "Eukaryotic RNA") + 
  my_plot_themes
PieChart_euk


ggsave(PieChart_euk,
       file = "Eukaryotic_RNA_PieChart.pdf",
       path = "Pie_Figures",
       width = 6, height = 4, units = "in")



###########################################################
##################### BACTERIAL PIE #######################


# Define custom colors for RNA_Type
bact_rna_colors <- c("bacterial rRNA" = "#FBE0B1", 
                "bacterial tRNA/ncRNA/etc" = "#FBB75B", 
                "bacterial mRNA" = "#D95319")

bact_data <- data.frame(
  RNA_Type = factor(c("bacterial mRNA", "bacterial rRNA", "bacterial tRNA/ncRNA/etc"), 
                    levels = c("bacterial mRNA", "bacterial rRNA", "bacterial tRNA/ncRNA/etc")),
  Percent = c(5, 85, 10)
)

PieChart_bact <- bact_data %>% 
  arrange(RNA_Type, desc(Percent)) %>%  # Need this so the numbers go to the correct slices
  group_by(RNA_Type) %>% # Need this or it won't be a round pie! 
  ungroup() %>%
  ggplot(aes(x = "", y = Percent, fill = RNA_Type)) +
  geom_bar(width = 1, stat = "identity", color = "black") + 
  coord_polar(theta = "y", start = 0) + 
  scale_fill_manual(values = bact_rna_colors) + 
  geom_text(aes(label = paste0(Percent, "%")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black") + 
  labs(title = "Bacterial RNA") + 
  my_plot_themes
PieChart_bact

ggsave(PieChart_bact,
       file = "Bacterial_RNA_PieChart.pdf",
       path = "Pie_Figures",
       width = 6, height = 4, units = "in")





