library(shiny)
library(tidyverse)
library(ggplot2)
library(minpack.lm)
library(ggpmisc)


source('../functions/pre_processing.R')
source('../functions/data_summary.R')
source('../functions/characterization.R')

# UI
ui <- fluidPage(
  titlePanel("Graphic user interface v1.1"),
  
  tabsetPanel(
    tabPanel("Main",
             tags$style(HTML("
        .tab-content {
          display: grid;
          grid-template-columns: 1fr;
          grid-template-rows: auto auto auto auto;
          gap: 10px;
          height: 100vh;
        }
        .centered-container {
          grid-row: 2;
          text-align: center;
        }
        .image-container {
          grid-row: 3;
          text-align: center;
        }
        .coded-by-container {
          grid-row: 4;
          text-align: center;
        }
      ")),
             div(class = "centered-container",
                 h2("Biosynthesis allocation model")
             ),
             div(class = "centered-container",
                 tags$p("Biosynthesis allocation model (BAM) is a graphic user interface based on the work of Nicolas Arias, Dahlin Zevallos-Aliaga, Tom Peeters and Daniel Guerra, available at",
                        tags$a("bioRxiv", href = "https://www.biorxiv.org/content/10.1101/2023.12.01.569588v1"))
             ),
             div(class = "image-container",
                 tags$img(src = "flux_allocation.png", alt = "Illustration of flux allocation", height = "300px")
             ),
             div(class = "centered-container",
                 div(style = "display: inline-block;",
                     fileInput("file", "Upload data")
                 )
             ),
             div(class = "coded-by-container",
                 p("Coded by Frank Britto Bisso")
             )
    ),
    tabPanel("Data summary",
             sidebarPanel(
               radioButtons("visualization", "Visualization", choices = c("Face wrap", "Single"), selected = "Face wrap"),
               selectInput("model_params", "Model parameters", choices = c("Optical density", "Fluorescence", "Growth rate", "Production rate", "Heterologous fraction")),
               fluidRow(
                 column(6, numericInput("start", "Start", value = 0, width = "100%")),
                 column(6, numericInput("end", "End", value = 1, width = "100%"))
               )
             ),
             mainPanel(
               plotOutput("raw_data_plot")
             )
    ),
    tabPanel("Characterization",
             sidebarPanel(
               selectInput("method_selection", "Method selection", choices = c("Single-point method", "Slope-based method")),
               selectInput("regression_method", "Regression method", choices = c("nls (nonlinear least squares)")),
               numericInput("H", "Maximum expression (H)", value = 8000),
               numericInput("ki", "Affinity for inducer (ki)", value = 50),
               numericInput("n", "Cooperativity (n)", value = 2),
               actionButton("fitting_button", "Fitting"),
               actionButton("restart_button", "Restart"),  # Add the "Restart" button
               br(),
               br(),
               selectInput("experimental_condition", "Experimental condition", choices = NULL),
               # Conditionally show t1, t2, and commit_changes_button only when "Slope-based method" is selected
               conditionalPanel(
                 condition = 'input.method_selection == "Slope-based method"',
                 fluidRow(
                   column(6, numericInput("t1", "t1", value = 0)),
                   column(6, numericInput("t2", "t2", value = 0)),
                   actionButton("commit_changes_button", "Commit changes")
                 )
               )
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Hill fitting",
                          mainPanel(
                            plotOutput("dose_plot")
                          )
                 ),
                 tabPanel("Quality control",
                          mainPanel(
                            fluidRow(
                              column(6, plotOutput("quality_control_plot1")),
                              conditionalPanel(
                                condition = 'input.method_selection == "Slope-based method"',
                                column(6, verbatimTextOutput("quality_control_table"))
                              ),
                              conditionalPanel(
                                condition = 'input.method_selection != "Slope-based method"',
                                column(6, plotOutput("quality_control_plot2"))
                              )
                            ),
                            fluidRow(
                              column(6, plotOutput("quality_control_plot3")),
                              column(6, plotOutput("quality_control_plot4"))
                            )
                          )
                 )
               )
             )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values
  uploaded_data <- reactiveVal(NULL)
  end_reactive <- reactiveVal(0)
  raw_data_plot <- reactiveVal(NULL)
  time_limits_reactive <- reactive(c(input$start, input$end))
  initial_guess_reactive <- reactive(list(H = input$H, ki = input$ki, n = input$n))
  file_uploaded <- reactiveVal(FALSE)
  experimental_condition_reactive <- reactiveVal(NULL)
  linear_detection_results_reactive <- reactiveVal(NULL)
  
  # Reactive values for tracking the "Fitting" button clicks
  fitting_button_single_point <- reactiveVal(FALSE)
  fitting_button_slope_based <- reactiveVal(FALSE)
  
  # Render the dose-response plot
  output$dose_plot <- renderPlot({
    req(file_uploaded())
    
    if (input$method_selection == "Single-point method") {
      dose_response_result <- dose_response_visualization(
        name = uploaded_data(),
        time_limits = time_limits_reactive(),
        initial_guess = initial_guess_reactive(),
        flag = fitting_button_single_point()
      )
      # Return the ggplot object from the list
      dose_response_result[[1]]
    } else if (input$method_selection == "Slope-based method") {
      slope_based_result <- slope_based_visualization(
        linear_df = linear_detection_results_reactive(),
        initial_guess = list(H = input$H, ki = input$ki, n = input$n),
        flag = fitting_button_slope_based()
      )
      # Return the ggplot object from the list
      slope_based_result[[1]]
    }
  })
  
  # Update the dataframe when a file is uploaded
  observeEvent(input$file, {
    uploaded_data(read.csv(input$file$datapath))
    file_uploaded(TRUE)
    
    num_unique_conditions <- length(unique(uploaded_data()$Condition))
    time_values <- unique(uploaded_data()$time)
    
    corresponding_time <- time_values[ceiling(nrow(uploaded_data()) / num_unique_conditions * 100) / 100]
    corresponding_time <- round(corresponding_time, 2)
    
    end_reactive(corresponding_time)
    updateNumericInput(session, "end", value = end_reactive())
    
    raw_data_plot(visualize_raw_data(
      data = uploaded_data(),
      x_var = "time",
      y_var = switch(input$model_params,
                     "Optical density" = "od",
                     "Fluorescence" = "flu",
                     "Heterologous fraction" = "phi",
                     "Production rate" = "production_rate",
                     "Growth rate" = "growth_rate"
      ),
      visualization = if (input$visualization == "Face wrap") "multiple" else "single",
      start = input$start,
      end = input$end
    ))
    
    # Update choices for the Experimental condition selectInput
    updateSelectInput(session, "experimental_condition", choices = unique(uploaded_data()$Condition))
    
    # Generate linear detection results once when data is uploaded
    linear_detection_results <- generate_linear_detection_results(uploaded_data())
    
    # Save the results in a reactiveVal for later use
    linear_detection_results_reactive(linear_detection_results)
  })
  
  # Update the plot when "Start" or "End" changes
  observeEvent(c(input$start, input$end), {
    req(file_uploaded())
    raw_data_plot(visualize_raw_data(
      data = uploaded_data(),
      x_var = "time",
      y_var = switch(input$model_params,
                     "Optical density" = "od",
                     "Fluorescence" = "flu",
                     "Heterologous fraction" = "phi",
                     "Production rate" = "production_rate",
                     "Growth rate" = "growth_rate"
      ),
      visualization = if (input$visualization == "Face wrap") "multiple" else "single",
      start = input$start,
      end = input$end
    ))
  })
  
  # Update the plot when "Visualization" changes
  observeEvent(input$visualization, {
    req(file_uploaded())
    raw_data_plot(visualize_raw_data(
      data = uploaded_data(),
      x_var = "time",
      y_var = switch(input$model_params,
                     "Optical density" = "od",
                     "Fluorescence" = "flu",
                     "Heterologous fraction" = "phi",
                     "Production rate" = "production_rate",
                     "Growth rate" = "growth_rate"
      ),
      visualization = if (input$visualization == "Face wrap") "multiple" else "single",
      start = input$start,
      end = input$end
    ))
  })
  
  observeEvent(input$model_params, {
    req(file_uploaded())
    y_var <- switch(input$model_params,
                    "Optical density" = "od",
                    "Fluorescence" = "flu",
                    "Heterologous fraction" = "phi",
                    "Production rate" = "production_rate",
                    "Growth rate" = "growth_rate"
    )
    
    visualization_type <- if (input$visualization == "Face wrap") {
      "multiple"
    } else {
      "single"
    }
    
    raw_data_plot(visualize_raw_data(
      data = uploaded_data(),
      x_var = "time",
      y_var = y_var,
      visualization = visualization_type,
      start = input$start,
      end = input$end
    ))
  })
  
  # Render the ggplot in the main panel
  output$raw_data_plot <- renderPlot({
    req(file_uploaded())
    raw_data_plot()
  })
  
  # Observe "Fitting" button clicks
  observeEvent(input$fitting_button, {
    req(file_uploaded())
    
    if (input$method_selection == "Single-point method") {
      fitting_button_single_point(TRUE)
    } else if (input$method_selection == "Slope-based method") {
      fitting_button_slope_based(TRUE)
    }
  })
  
  # Observe "Restart" button clicks
  observeEvent(input$restart_button, {
    # Reset only the reactive values associated with the "Fitting" button
    fitting_button_single_point(FALSE)
    fitting_button_slope_based(FALSE)
  })
  
  observeEvent(input$experimental_condition, {
    req(file_uploaded())
    
    # Check if linear_detection_results_reactive is not NULL
    if (!is.null(linear_detection_results_reactive())) {
      
      # Find the corresponding row in linear_detection_results
      condition_row <- which(linear_detection_results_reactive()$Condition == input$experimental_condition)
      
      # Check if a matching row is found
      if (length(condition_row) > 0) {
        # Update t1 and t2 based on the corresponding values in linear_detection_results
        updateNumericInput(session, "t1", value = linear_detection_results_reactive()$Minimum[condition_row])
        updateNumericInput(session, "t2", value = linear_detection_results_reactive()$Maximum[condition_row])
      } else {
        # Handle the case where no matching row is found (optional)
        # You might want to set default values or take some other action
        # For example:
        updateNumericInput(session, "t1", value = 0)
        updateNumericInput(session, "t2", value = 0)
      }
    }      
  })
  
  # Observe changes in "experimental_condition" and "method_selection"
  observeEvent(c(input$experimental_condition, input$method_selection, input$t1, input$t2), {
    req(file_uploaded())
    
    # Check if linear_detection_results_reactive is not NULL
    if (!is.null(linear_detection_results_reactive())) {
      
      # Additional logic based on "Single-point selection"
      if (input$method_selection == "Single-point method") {
        # Call the "accumulation" function
        accumulation_plots <- accumulation(
          name = uploaded_data(),
          inductor = input$experimental_condition
        )
        
        # Display the resulting ggplot2 objects in the main panel of the "Quality control" sub-panel
        output$quality_control_plot1 <- renderPlot({
          accumulation_plots[[1]]
        }, width = 350, height = 400)
        output$quality_control_plot2 <- renderPlot({
          # Display the second element of the list in the "Hill fitting" tab panel
          if (fitting_button_slope_based() || fitting_button_single_point()) {
            accumulation_plots[[2]]
          }
        }, width = 350, height = 400)
        
        # Hide the other plots
        output$quality_control_plot3 <- NULL
        output$quality_control_plot4 <- NULL
        output$quality_control_table <- NULL
        
      }
      
      # Additional logic based on "Method selection"
      else if (input$method_selection == "Slope-based method") {
        
        # Call the "hook_plot" function
        hook_plot_output <- hook_plot(
          name = uploaded_data(),
          time_limits = list(start = input$start, end = input$end)
        )
        
        # Call the "visualize_manual_selection" function
        manual_selection_plots <- visualize_manual_selection(
          name = uploaded_data(),
          t1 = input$t1,
          t2 = input$t2,
          condition = input$experimental_condition
        )
        
        # Render the plots in the Shiny app
        output$quality_control_plot1 <- renderPlot({
          hook_plot_output
        }, width = 300, height = 300)
        
        output$quality_control_plot2 <- renderPlot({
          # Display the second element of the list in the "Hill fitting" tab panel
          if (fitting_button_slope_based() || fitting_button_single_point()) {
            manual_selection_plots[[2]]
          }
        }, width = 300, height = 300)
        
        output$quality_control_table <- renderPrint({
          linear_detection_results_reactive()
        })
        
        output$quality_control_plot3 <- renderPlot({
          manual_selection_plots[[1]]
        }, width = 300, height = 300)
        
        output$quality_control_plot4 <- renderPlot({
          manual_selection_plots[[2]]
        }, width = 300, height = 300)
      }
    }
  })
  
  # Observe changes when the "Commit changes" button is clicked
  observeEvent(input$commit_changes_button, {
    req(file_uploaded())
    
    # Check if linear_detection_results_reactive is not NULL
    if (!is.null(linear_detection_results_reactive())) {
      
      # Get the current values of "t1" and "t2"
      t1_value <- input$t1
      t2_value <- input$t2
      
      # Get the current value of "Experimental condition"
      current_condition <- input$experimental_condition
      
      # Find the corresponding row in linear_detection_results
      condition_row <- which(linear_detection_results_reactive()$Condition == current_condition)
      
      # Check if a matching row is found
      if (length(condition_row) > 0) {
        # Update "Minimum" and "Maximum" based on the current values of "t1" and "t2"
        linear_detection_results_data <- linear_detection_results_reactive()
        linear_detection_results_data$Minimum[condition_row] <- t1_value
        linear_detection_results_data$Maximum[condition_row] <- t2_value
        linear_detection_results_reactive(linear_detection_results_data)
      } else {
        # Handle the case where no matching row is found (optional)
        # You might want to set default values or take some other action
        # For example:
        cat("Warning: No matching row found for the current condition.")
      }
    }
  })
}


# Run the Shiny app
shinyApp(ui, server)



