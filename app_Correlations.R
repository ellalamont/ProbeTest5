# 3/3/25
# Interactive plot of the corrlation between the average TPM scaled values for THP1 cells spiked with 1e6 cells oF H37Ra.
# Doing this so I can look at where the gene sets are in the correlation!

library(shiny)
library(gmodels)

source("Import_data.R") 

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12), 
        plot.subtitle = element_text(size=12), 
        plot.margin = margin(10, 10, 10, 20))

# Process the data a little more (not done in import data)
ProbeTest5_tpm_2 <- ProbeTest5_tpm
# Adjust the names so they are slightly shorter
names(ProbeTest5_tpm_2) <- gsub(x = names(ProbeTest5_tpm_2), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it
# add rownames to the tpm and metadata dataframes
rownames(ProbeTest5_tpm_2) <- ProbeTest5_tpm_2[,1] # add the rownames
ProbeTest5_tpm_2 <- ProbeTest5_tpm_2[,-1] # Remove the old column of rownames
# Need to make sure there is a Gene column (gets lost)
ProbeTest5_tpm_2$Gene <- rownames(ProbeTest5_tpm_2)

# Log10 transform the data
ProbeTest5_tpm_Log10 <- ProbeTest5_tpm_2 %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values


# Take the average of the spiked capture and the averages of the spiked not capture and see if that is any better

THP1_AVERAGE_tpm_Log10 <- ProbeTest5_tpm_Log10 %>% select(Gene, THP1_1e6_1a, THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_2b, THP1_1e6_3a, THP1_1e6_3b) %>%
  mutate(
    CAPTURED_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    NOTCaptured_THP1Spiked1e6 = rowMeans(select(., c(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b)), na.rm = TRUE),
  )

Sample1 <- "CAPTURED_THP1Spiked1e6" # Captured
Sample2 <- "NOTCaptured_THP1Spiked1e6" # Not Captured
ScatterCorr <- THP1_AVERAGE_tpm_Log10 %>% 
  # filter(NOTCaptured_THP1Spiked1e6 != 0) %>% # remove all the zeros in the not captured sample and increases the correlation
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation",
       x = paste0("Log10(TPM+1) CAPTURED samples averaged"), y = paste0("Log10(TPM+1) NOT captured samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr






# Define UI ----
ui <- fluidPage(
  titlePanel("THP1 spiked with 1e6cells of H37Ra"),
  
  fluidRow(
    
    column(width = 2.5,
           fluidRow(
             column(width = 3,
                    textInput("my_GeneID", 
                              label = "Gene ID (Will label gene orange)",
                              value = "Rv...")
             )
           ),
           selectInput("my_GeneSet",
                       label = "Gene Set (Will label genes yellow)",
                       choices = names(allGeneSets))
    ),
    
    column(width = 5,
           plotlyOutput("Correlation_plot",
                        width = 1000, height = 600),
    ),
  )
  
)


# Define server logic ----
server <- function(input, output) {
  
  # Correlation Plot
  
  output$Correlation_plot <- renderPlotly({
    
    # Add data for labelling a single gene
    single_gene <- THP1_AVERAGE_tpm_Log10 %>% 
      filter(Gene == input$my_GeneID)
    
    # Add data for labelling a gene set
    gene_set <- THP1_AVERAGE_tpm_Log10 %>%
      filter(Gene %in% allGeneSets[[input$my_GeneSet]])
    
    ScatterCorr <- THP1_AVERAGE_tpm_Log10 %>% 
      ggplot(aes(x = CAPTURED_THP1Spiked1e6, y = NOTCaptured_THP1Spiked1e6)) + 
      geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
      
      # Add a differently colored point
      geom_point(data = gene_set, color = "yellow",  aes(label = Gene)) + 
      geom_point(data = single_gene, color = "orange", aes(label = Gene)) + 
      
      geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
      labs(title = paste0("THP1 spiked w/1e6 H37Rv Samples averaged"),
           subtitle = "Pearson correlation",
           x = paste0("Log10(TPM+1) CAPTURED samples averaged"), y = paste0("Log10(TPM+1) NOT captured samples averaged")) + 
      stat_cor(method="pearson") + # add a correlation to the plot
      my_plot_themes
    ScatterCorr
    
    
    
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)