# Probe Prep testing 5
# Probes made by Jessica, run included THP1 test samples and sputum
# See Test_Sputum_Jan2025_RNALibraryPrep
# 2/15/25

# 5/9/25: Updated TPM and RPKM files to the ones that had KEEPINTERGNICS = F
# Also removing all the scaled stuff because not using scaled anymore
# All the old figure folders have been renamed to include _old.KIT (KEEP INTERGENICS T)


################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(scales) # To add commas to the y axis... scale_y_continuous(label=comma, )

################################################
################### Colors #####################

cbPalette_1 <- c("#999999", "#E69F00") # Gold and Grey
cbPalette_1.5 <- c("#E69F00", "#999999") # Gold and Grey
cbPalette_2 <- c( "#0072B2", "#999999") # Blue and Grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 <-  c("#bfbfbf", "#56B4E9")
cbPalette3 <-  c("#bfbfbf", "#E69F00")
cbPalette4 <- c("#56B4E9", "#009E73", "#F0E442")
cbPalette5 <- c("#009E73", "#FF7F00")
c25 <- c(
  "dodgerblue2", "#E31A1C", "green4",
  "#6A3D9A","#FF7F00","black", "gold1",
  "skyblue2", "#FB9A99","palegreen2","#CAB2D6",
  "#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown"
)
c12 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "palegreen2", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4") 
c16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black","gold1", "#FB9A99", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4") 
c11 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "palegreen2", "gold1", "maroon", "orchid1", "darkturquoise", "darkorange4", "gray70")
c7 <- c("gray70", "#E0F7FA", "#B2EBF2", "#81D4FA", "#03A9F4","#0288D1", "#01579B")
c8 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "palegreen2", "gray70", "maroon", "black")
Orange_Gradient <- c("#FFD4A3", "#FFA64D", "#FF7F00", "#E66900", "#993D00")
c14 <- c("dodgerblue2", "#E31A1C", "green4", "#FF7F00", "black","gold1", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4") 
c3 <- c("#56B4E9", "#E66900", "#009E73")

# Stop scientific notation
options(scipen = 999) 
# options(scipen = 0) # To revert back to default

###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############
# Importing the ProbeTests 3 and 4 and 5 to get all the sputum samples I have done

# ProbeTest5_pipeSummary <- read.csv("ProbeTest5_Pipeline.Summary.Details.csv")
# This has been edited to include more metadata!
ProbeTest5_pipeSummary <- read.csv("ProbeTest5_Pipeline.Summary.Details_moreTrim.csv") # This has the 3' end trimmed 40bp to increase the number of reads aligning
ProbeTest4_pipeSummary <- read.csv("ProbeTest4_Pipeline.Summary.Details.csv")
ProbeTest3_pipeSummary <- read.csv("ProbeTest3_Pipeline.Summary.Details.csv")

# Merge the 3 documents
All_pipeSummary <- merge(ProbeTest5_pipeSummary, ProbeTest4_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, ProbeTest3_pipeSummary, all = T)

All_pipeSummary$X <- NULL
All_pipeSummary$X.1 <- NULL

All_pipeSummary$Hyb_Time <- as.character(All_pipeSummary$Hyb_Time)
ordered_Hyb_Time <- c("4", "16")
All_pipeSummary$Hyb_Time <- factor(All_pipeSummary$Hyb_Time, levels = ordered_Hyb_Time)

All_pipeSummary$Week <- as.character(All_pipeSummary$Week)
ordered_Week <- c("0", "2", "4")
All_pipeSummary$Week <- factor(All_pipeSummary$Week, levels = ordered_Week)

All_pipeSummary$EukrRNADep <- as.character(All_pipeSummary$EukrRNADep)
ordered_EukrRNADep <- c("MtbrRNA", "DualrRNA")
All_pipeSummary$EukrRNADep <- factor(All_pipeSummary$EukrRNADep, levels = ordered_EukrRNADep)

# Remove the undetermined
All_pipeSummary <- All_pipeSummary %>% filter(SampleID != "Undetermined_S0")

# Remove the marmoset and the high low THP1 samples
All_pipeSummary <- All_pipeSummary %>% filter(!Sample_Type %in% c("Marmoset", "High_Low_THP1"))

All_pipeSummary$SampleID <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

All_pipeSummary <- All_pipeSummary %>% mutate(Sputum_Number = str_extract(SampleID, "S_[0-9]+"))


###########################################################
####### PIPE SUMMARY: EXTRACT JUST THE SPUTUM SAMPLES #####

# NOT including the DualrRNA!!

AllSputum_pipeSummary <- All_pipeSummary %>% filter(Sample_Type == "Sputum")

# Just keep unique sputum samples, the best sequencing results
# SampleID to Keep:
# Unique Sputum: 
# W0 samples: "S_250754", "S_355466", "S_503557" 
# W2 samples: "S_349941_Probe_3D_25", "S_503937", "S_575533_MtbrRNA", "S_577208"
# W4 samples: "S_351946_Probe_4A_100", "S_687338_Probe_4A_100"

Unique_Sputum <- c("S_250754", "S_355466", "S_503557", "S_349941_Probe_3D_25", "S_503937", "S_575533_MtbrRNA", "S_577208", "S_351946_Probe_4A_100", "S_687338_Probe_4A_100")

UniqueSputum_pipeSummary <- AllSputum_pipeSummary %>% filter(SampleID %in% Unique_Sputum)

###############################################################
####### PIPE SUMMARY: EXTRACT LIMIT OF DETECTION SAMPLES ######

LimitofDetect_pipeSummary <- All_pipeSummary %>% 
  filter(Sample_Type == "THP1") %>% 
  filter(Ra_cells != "none") %>% 
  filter(Probe != "None") %>%
  mutate(Ra_cells2 = case_when(
    Ra_cells == "one_e_2" ~ "1e2",
    Ra_cells == "one_e_3" ~ "1e3",
    Ra_cells == "one_e_4" ~ "1e4",
    Ra_cells == "one_e_5" ~ "1e5",
    Ra_cells == "one_e_6" ~ "1e6",
    Ra_cells == "one_e_8" ~ "1e8",
    TRUE ~ "default_value"  # this is optional for values that don't meet any condition
  ))
  
# Should maybe go back and remove all the 1e6 cells that are not specifically for the limit of detection, but not doing that right now.
# Or what might be easier is just include the samples from ProbeTest5, where I have at least 3 replicates each! 



###############################################################
######## PIPE SUMMARY: EXTRACT CAPTURED VS NOT SAMPLES ########

# Just get my 6 samples that are matchine
CapturedVsNot_pipeSummary <- All_pipeSummary %>%
  filter(Sample_Type == "THP1") %>% 
  filter(Ra_cells == "one_e_6") %>% 
  filter(Run == "ProbeTest5") %>%
  filter(SampleID != "THP1_1e6_4") %>%
  mutate(Replicates = c("R1","R1","R2","R2","R3", "R3"))
  
ordered_Probe <- c("None", "JA2")
CapturedVsNot_pipeSummary$Probe <- factor(CapturedVsNot_pipeSummary$Probe, levels = ordered_Probe)

###########################################################
############ IMPORT AND PROCESS ALL TPM VALUES ############

# # ProbeTest5_tpm <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.SCALED.TPM.csv")
# ProbeTest5_tpm <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.SCALED.TPM_moreTrim.csv") # This has the 3' end trimmed 40bp to increase the number of reads aligning
# ProbeTest4_tpm <- read.csv("ProbeTest4_Mtb.Expression.Gene.Data.SCALED.TPM.csv")
# ProbeTest3_tpm <- read.csv("ProbeTest3_Mtb.Expression.Gene.Data.SCALED.TPM.csv")
# 
# # Need to remove the undetermined which all share names
# ProbeTest5_tpm$Undetermined_S0 <- NULL
# ProbeTest4_tpm$Undetermined_S0 <- NULL
# ProbeTest3_tpm$Undetermined_S0 <- NULL
# 
# # Merge the 3 documents
# All_tpm <- merge(ProbeTest5_tpm, ProbeTest4_tpm, all = T)
# All_tpm <- merge(All_tpm, ProbeTest3_tpm, all = T)
# 
# # Adjust the names so they are slightly shorter
# names(All_tpm) <- gsub(x = names(All_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)
# 
# # Grab the metadata I added to my_pipeSummary
my_metadata <- All_pipeSummary %>% select(1, 13:ncol(All_pipeSummary))
# 
# # add rownames to the tpm and metadata dataframes
# rownames(All_tpm) <- All_tpm[,1] # add the rownames
# All_tpm <- All_tpm[,-1] # Remove the old column of rownames
# rownames(my_metadata) <- my_metadata[,1] # add the rownames
# my_metadata <- my_metadata[,-1] # Remove the old column of rownames

# Error in `.rowNamesDF<-`(x, value = value) : 
#   duplicate 'row.names' are not allowed
# In addition: Warning message:
#   non-unique value when setting 'row.names': ‘THP1_1e6_4’ 


###########################################################
########## TPM: EXTRACT JUST THE SPUTUM SAMPLES ###########

# # Combine the sputum samples only 
# UniqueSputum_tpm <- All_tpm %>% select(all_of(Unique_Sputum))
# # write.csv(AllSputum_tpm, "AllSputum_tpm.csv")
# 
# # Grab AllSputum metadata
UniqueSputum_metadata <- UniqueSputum_pipeSummary # %>% select(2, 14:30)

# Adjust the metadata names so they are the same
rownames(UniqueSputum_metadata) <- UniqueSputum_metadata[,1] # add the rownames

# Get the names of the samples with >1M reads
Unique_Sputum_1Mreads <- UniqueSputum_metadata %>% filter(N_Genomic > 1000000) %>% pull(SampleID)



###########################################################
###### IMPORT AND PROCESS ALL NOT scaled TPM VALUES #######

# ProbeTest5_tpm <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.SCALED.TPM.csv")
ProbeTest5_tpm_NOTscaled <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.TPM_moreTrim.csv") # This has the 3' end trimmed 40bp to increase the number of reads aligning
ProbeTest4_tpm_NOTscaled <- read.csv("ProbeTest4_Mtb.Expression.Gene.Data.TPM.csv")
ProbeTest3_tpm_NOTscaled <- read.csv("ProbeTest3_Mtb.Expression.Gene.Data.TPM.csv")

# Need to remove the undetermined which all share names
ProbeTest5_tpm_NOTscaled$Undetermined_S0 <- NULL
ProbeTest4_tpm_NOTscaled$Undetermined_S0 <- NULL
ProbeTest3_tpm_NOTscaled$Undetermined_S0 <- NULL

# Merge the 3 documents
All_tpm_NOTscaled <- full_join(ProbeTest5_tpm_NOTscaled, ProbeTest4_tpm_NOTscaled, by = "X")
All_tpm_NOTscaled <- full_join(All_tpm_NOTscaled, ProbeTest3_tpm_NOTscaled, by = "X")
All_tpm_NOTscaled <- All_tpm_NOTscaled %>% rename("Gene" = "X")

# Adjust the names so they are slightly shorter
names(All_tpm_NOTscaled) <- gsub(x = names(All_tpm_NOTscaled), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# add rownames to the tpm and metadata dataframes
rownames(All_tpm_NOTscaled) <- All_tpm_NOTscaled[,1] # add the rownames


############################################################
###### NOT SCALED TPM: EXTRACT JUST THE SPUTUM SAMPLES #####

# Combine the sputum samples only 
UniqueSputum_tpm_NOTscaled <- All_tpm_NOTscaled %>% select(all_of(Unique_Sputum))
# write.csv(AllSputum_tpm, "AllSputum_tpm.csv")

###########################################################
############ TPM: ADJUST PROBETEST5 NOT SCALED ############

# Adjust the names so they are slightly shorter
names(ProbeTest5_tpm_NOTscaled) <- gsub(x = names(ProbeTest5_tpm_NOTscaled), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)
rownames(ProbeTest5_tpm_NOTscaled) <- ProbeTest5_tpm_NOTscaled[,1] # add the rownames
ProbeTest5_tpm_NOTscaled <- ProbeTest5_tpm_NOTscaled %>% rename(Gene = X)

###########################################################
####### TPM: EXTRACT CAPTURED VS NOT THP1 SAMPLES #########

# Want the captured samples to be scaled and the not captured samples to be not scaled
# Captured = "THP1_1e6_1a", "THP1_1e6_2b", "THP1_1e6_3a" 
# Not Captured = "THP1_1e6_1b", "THP1_1e6_2a", "THP1_1e6_3b" 

# # Adjust the names so they are slightly shorter
# names(ProbeTest5_tpm) <- gsub(x = names(ProbeTest5_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it
# # add rownames to the tpm and metadata dataframes
# rownames(ProbeTest5_tpm) <- ProbeTest5_tpm[,1] # add the rownames
# ProbeTest5_tpm <- ProbeTest5_tpm[,-1] # Remove the old column of rownames
# # Need to make sure there is a Gene column (gets lost)
# ProbeTest5_tpm$Gene <- rownames(ProbeTest5_tpm)
# 
# THP1Spike_tpm_CorrectScales <- inner_join(
#   ProbeTest5_tpm %>% select(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a, Gene),
#   ProbeTest5_tpm_NOTscaled %>% select(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b, Gene),
#   by = "Gene"
# )

###########################################################
############# IMPORT AND PROCESS ALL READS_M ##############
# 5/1/25

# ProbeTest5_tpm <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.SCALED.TPM.csv")
ProbeTest5_RawReads <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.ReadsM_moreTrim.csv") # This has the 3' end trimmed 40bp to increase the number of reads aligning
ProbeTest4_RawReads <- read.csv("ProbeTest4_Mtb.Expression.Gene.Data.ReadsM.csv")
# Not doing ProbeTest3 because I didn't pull the not scaled from the lenovo and it doesn't contain the sputum samples >1M reads

# Need to remove the undetermined which all share names
ProbeTest5_RawReads$Undetermined_S0 <- NULL
ProbeTest4_RawReads$Undetermined_S0 <- NULL

# Merge the 3 documents
All_RawReads <- full_join(ProbeTest5_RawReads, ProbeTest4_RawReads, by = "X")
All_RawReads <- All_RawReads %>% rename("Gene" = "X")

# Adjust the names so they are slightly shorter
names(All_RawReads) <- gsub(x = names(All_RawReads), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# add rownames to the tpm and metadata dataframes
rownames(All_RawReads) <- All_RawReads[,1] # add the rownames


# Get just the sputum we are interested in and the not-captured broth
SputumBroth_RawReads <- All_RawReads %>% select(all_of(c(Unique_Sputum_1Mreads, "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6", "Gene")))




###########################################################
############ RPKM and R_M: IMPORT PROBETEST5 ##############
# Want to see if the broth correlations look different for RPKM or Reads_M values, using scaled for all.

# ProbeTest5_RPKM <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.SCALED.RPKM_moreTrim.csv")
# ProbeTest5_RPKM$Undetermined_S0 <- NULL
# # Adjust the names so they are slightly shorter
# names(ProbeTest5_RPKM) <- gsub(x = names(ProbeTest5_RPKM), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)
# rownames(ProbeTest5_RPKM) <- ProbeTest5_RPKM[,1] # add the rownames
# ProbeTest5_RPKM <- ProbeTest5_RPKM %>% rename(Gene = X)
# 
# Broth_RPKM <- ProbeTest5_RPKM %>% select(contains("Broth"))

# ProbeTest5_ReadsM <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.SCALED.ReadsM_moreTrim.csv")
# ProbeTest5_ReadsM$Undetermined_S0 <- NULL
# # Adjust the names so they are slightly shorter
# names(ProbeTest5_ReadsM) <- gsub(x = names(ProbeTest5_ReadsM), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)
# rownames(ProbeTest5_ReadsM) <- ProbeTest5_ReadsM[,1] # add the rownames
# ProbeTest5_ReadsM <- ProbeTest5_ReadsM %>% rename(Gene = X)
# 
# Broth_ReadsM <- ProbeTest5_ReadsM %>% select(contains("Broth"))

###########################################################
########### TPM: EXTRACT JUST THE BROTH SAMPLES ###########

# # Keep the broth samples only
# Broth_tpm <- All_tpm %>% select(contains("Broth"))
# 
# # Grab AllSputum metadata
# Broth_metadata <- All_pipeSummary %>% 
#   filter(grepl("Broth", SampleID)) %>%
#   select(where(~!all(is.na(.)))) %>% # Removing all the columns that are all NA
#   mutate(Probe = str_replace_all(Probe, c("JA2"="Captured", "None"="Not captured"))) # Changing what's in the Probe column
# rownames(Broth_metadata) <- Broth_metadata[,1] # add the rownames
# 
# # Do the same for the NOTscaled tpm
# Broth_tpm_NOTscaled <- ProbeTest5_tpm_NOTscaled %>% select(contains("Broth"))
# Broth_tpm_NOTscaled$Gene <- rownames(Broth_tpm_NOTscaled)
# 
# Broth_tpm$Gene <- rownames(Broth_tpm)
# 
# # Realized that I am using the not captured as scaled, but it doesn't need to be scaled! Make everything with the scaled captured and the not scaled not captured
# Broth_tpm_CorrectScales <- inner_join(
#   Broth_tpm %>% select(H37Ra_Broth_1, H37Ra_Broth_2, H37Ra_Broth_3, Gene),
#   Broth_tpm_NOTscaled %>% select(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6, Gene),
#   by = "Gene"
# )

###########################################################
############# SPUTUM ALIGN TO MULTIPLE GENOMES ############

# This isn't actually any samples from ProbeTest5, work was done on lenovo under EL_UniqueSputum, but doing it here because this is where I'm working right now

Diff_alignments <- read.csv2("EL_DifferentAlignments_Summary_byhand.csv", sep = ",")

# order the Genome_Alignment
ordered_Genome <- c("H37Rv", "Clin1", "CG24", "Mada116")
Diff_alignments$Genome_Alignment <- factor(Diff_alignments$Genome_Alignment, levels = ordered_Genome)


###########################################################
#################### IMPORT GENE SETS #####################


# load("MTb.KEGG.Pathways.rda")
# names(allGeneSets) <- gsub("<.*", "", names(allGeneSets))

load("MTb.Tuberculist.FunctionalGroups.rda")

# load("MTb.Biocyc.Pathways.rda")
