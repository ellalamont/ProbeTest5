# Compare summary data for 2 samples that were aligned to different clinical genomes
# E. Lamont
# 2/18/25

source("Import_data.R") # to get Diff_alignments

# From Bob
# 1) I made target genomes for 3 clinical isolates:  the first one you selected I called "Clin1", and then both "CG24" and "Mala116" from your later suggestions.  In the folder called "EL_UniqueSputum", I ran 2 of your samples against all 4 genomes {3 isolates and H37} using identical alignment settings (no ribo clearing, keep only real genes, K=2 for reads hitting multiple locations).  So you will find 4 different "results" folders with the isolate name explicit in the folder name.  Under each is the transcriptome you get for each sample, and look also at the "summary/xxx.pipelineSummary.txt" files for how many of the reads aligned to each isolate genome.  To my eye all 4 are very similar.  I am hoping you can do the "how many genes got at least 10 reads?" steps on your own?...  this should start to address your questions about "what might we be missing..."

# Note: The N_Genomic is including all the rRNA genes as well, this is why it is a larger number than before!

# Need to know how many genes are in each of these clinical genomes! 
# Made TPM files for each and looked at number of rows
# Clin1: 4155
# CG24: 4225
# Mada116: 4146

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=14),
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

facet_themes <- theme(strip.background=element_rect(fill="#ededed", size = 0.9),
                      strip.text = element_text(size = 10))


###########################################################
#################### >=10 READS ALIGNING ##################

AtLeast10Reads_v1 <- Diff_alignments %>%
  ggplot(aes(x = SampleID, y = AtLeast.10.Reads, fill = Genome_Alignment)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values=c("#FF5733", "#87CEFA", "#1E90FF", "#6495ED")) +
  scale_y_continuous(limits = c(0,4499), breaks = seq(0, 4500, 500)) + 
  geom_text(aes(label = AtLeast.10.Reads),
            position = position_dodge(width = 0.9),
            vjust = -0.5,  # adjust vertical position above the bar
            size = 3) +
  labs(title = "Sputum samples aligned to different clinical genomes",
       x = "Sputum sample",
       y = "# genes with at least 10 reads") +
  my_plot_themes
AtLeast10Reads_v1
ggsave(AtLeast10Reads_v1,
       file = "DiffAlignments_10Genes.pdf",
       path = "Diff_alignments_Figures",
       width = 7, height = 4, units = "in")

# Make it stacked given the total reads I think there are in each gene
AtLeast10Reads_v2 <- Diff_alignments %>% 
  mutate(AtLeast.10.Reads_Percent = case_when(
    Genome_Alignment == "H37Rv" ~ AtLeast.10.Reads/4499 * 100,
    Genome_Alignment == "Clin1" ~ AtLeast.10.Reads/4155 * 100,
    Genome_Alignment == "CG24" ~ AtLeast.10.Reads/4225 * 100,
    Genome_Alignment == "Mada116" ~ AtLeast.10.Reads/4146 * 100
  )) %>% 
  mutate(LessThan.10.Reads_Percent = case_when(
    Genome_Alignment == "H37Rv" ~ 100 - AtLeast.10.Reads_Percent,
    Genome_Alignment == "Clin1" ~ 100 - AtLeast.10.Reads_Percent,
    Genome_Alignment == "CG24" ~ 100 - AtLeast.10.Reads_Percent,
    Genome_Alignment == "Mada116" ~ 100 - AtLeast.10.Reads_Percent
  )) %>% 
  pivot_longer(cols = c("AtLeast.10.Reads_Percent", "LessThan.10.Reads_Percent"),
               names_to = "Read_Type", 
               values_to = "N_genes") %>% 
  mutate(Read_Type = factor(Read_Type, levels = c("LessThan.10.Reads_Percent", "AtLeast.10.Reads_Percent"))) %>% 
  ggplot(aes(x = Genome_Alignment, y = N_genes, fill = Read_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_text(aes(label = ifelse(Read_Type == "AtLeast.10.Reads_Percent", 
                               paste0(round(N_genes), "%"), 
                               "")),
            position = position_stack(vjust = 0.98),
            size = 4) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  facet_wrap(~ SampleID) +
  labs(title = "Sputum samples aligned to different clinical genomes",
       subtitle = "Got total number of genes from row number of TPM files",
       x = "Sputum sample",
       y = "% genes with at least 10 reads") +
  my_plot_themes + facet_themes
AtLeast10Reads_v2
ggsave(AtLeast10Reads_v2,
       file = "DiffAlignments_10Genes_Percent.pdf",
       path = "Diff_alignments_Figures",
       width = 8, height = 5, units = "in")


###########################################################
################## NUMBER READS ALIGNING ##################

Reads_aligning_v1 <- Diff_alignments %>%
  pivot_longer(cols = c("N_Genomic", "N_NoHit"),
               names_to = "Read_Type",
               values_to = "Read_Count") %>%
  mutate(Read_Type = factor(Read_Type, levels = c("N_NoHit", "N_Genomic"))) %>% 
  ggplot(aes(x = Genome_Alignment, y = Read_Count, fill = Read_Type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(limits = c(0,11000000), breaks = seq(0, 11000000, 1000000)) +
  scale_fill_manual(values=cbPalette_1) +
  facet_wrap(~ SampleID) +
  labs(title = "Sputum samples aligned to different clinical genomes",
       x = "Sputum sample",
       y = "# reads") +
  my_plot_themes + facet_themes
Reads_aligning_v1
ggsave(Reads_aligning_v1,
       file = "DiffAlignments_ReadsAligning.pdf",
       path = "Diff_alignments_Figures",
       width = 7, height = 4, units = "in")




