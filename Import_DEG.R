# Import Bob's DE Data
# E. Lamont
# 2/20/25

###########################################################
################## IMPORT BOBS DE DATA #################### 

# Naming these by hand which is a pain but not sure how else to do it to make sure I am keeping track of which one is being compared to which

AllBrothCaptured_ComparedTo_AllBrothNOTCaptured <- read.delim("JOINED_BobAverages/MTb.MetaResults.Broth_Captured_vs_NOTCaptured/Broth_captured.MTb.Meta.JOINED.txt")

W2_ComparedTo_W0 <- read.delim("JOINED_BobAverages/MTb.MetaResults.W2_vs_W0/W2.MTb.Meta.JOINED.txt")


###########################################################
################ MAKE A LIST OF ALL DFs ###################

list_dfs <- list(AllBrothCaptured_ComparedTo_AllBrothNOTCaptured, 
                 W2_ComparedTo_W0) 

# Make a list of all the names
df_names <- c("AllBrothCaptured_ComparedTo_AllBrothNOTCaptured", 
              "W2_ComparedTo_W0"
)

# Give the df list the correct df names
names(list_dfs) <- df_names




###########################################################
############### ADD COLUMNS OF DE VALUES ##################

# Make a new list to hold dataframes with extra columns
list_dfs_2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs)) {
  
  current_df <- list_dfs[[i]]
  current_df_name <- df_names[i]
  
  # Make the column pointing out which ones are differential expressed
  current_df$DE <- ifelse(current_df$LOG2FOLD < -2 & current_df$AVG_PVALUE < 0.05, "significant down",
                          ifelse(current_df$LOG2FOLD > 2 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE <- factor(current_df$DE, levels = ordered_DE)
  
  # Make the column with DE gene names for plotting on graph
  current_df$DE_labels <- ifelse(current_df$DE != "not significant", current_df$GENE_NAME, NA)
  
  list_dfs_2[[current_df_name]] <- current_df
  
}









