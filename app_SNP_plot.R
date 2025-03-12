library(shiny)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)

source("Import_Data_SNPs.R") # SNPs_df_4

###########################################################
################# MAKE USER INPUT LISTS ###################

sample_list <- SNPs_df_4$SampleID %>% as.factor() %>% levels

###########################################################

my_colors <- c("grey", "green", "orange", "blue", "red", "black", "purple")

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(angle = 45, size=10, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        plot.subtitle = element_text(size=10),
        plot.margin = margin(10, 10, 10, 20))


# Define UI ----
ui <- fluidPage(
  titlePanel("EL ProbeTest5 phoP (Rv0757) SNPs from captured samples"),
  
  fluidRow(
    
    column(width = 2,
           h3("Variables to change"),
           selectInput("my_sample",
                       label = "Choose sample",
                       choices = sample_list),
           radioButtons("Boxplot_y",
                        label = "Choose y axis",
                        choices = list("fill", "stack")),
           sliderInput("range",
                       label = "Choose position range",
                       min = 851608, max = 852351,
                       value = c(851608, 852351)),
    ),
    
    column(width = 5,
           plotlyOutput("stacked_bar",
                        width = 1300, height = 600),
    ),
  )
  
)

# Define server logic ----
server <- function(input, output) {
  
  output$stacked_bar <- renderPlotly({
    
    my_title <- paste0(input$my_sample)
    my_subtitle <- NULL
    my_x_axis <- "Position"
    
    my_y_axis <- switch(input$Boxplot_y,
                        "fill" = "Proportion of each base",
                        "stack" = "Raw Reads of each base")
    
    my_labels <- labs(title = my_title,
                      x = my_x_axis, 
                      y = my_y_axis,
                      subtitle = my_subtitle)
    
    my_barplot <- SNPs_df_4 %>%
      filter(SampleID == input$my_sample) %>%
      filter(between(POSITION, input$range[1], input$range[2])) %>% 
      ggplot(aes(x=POSITION, y=raw_read_number, fill=base_call, text = paste0("Reference Base: ", REF_BASE))) + 
      geom_bar(position=input$Boxplot_y, stat="identity") + # position = stack for raw reads, position = fill for proportion
      scale_fill_manual(values=my_colors) + 
      scale_x_continuous(breaks = seq(from = input$range[1], to = input$range[2], by = 2))
    
    final_plot <- my_barplot + my_plot_themes + my_labels
    
    final_plot_plotly <- ggplotly(final_plot)
    
    final_plot_plotly
    
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)