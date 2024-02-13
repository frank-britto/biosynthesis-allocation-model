# Load libraries
library(shiny)
library(shinyjs)
library(readxl)
library(dplyr)

# Import functions
source('../functions/pre_processing.R')

# Function to generate well names based on grid dimensions
generate_well_names <- function(rows, cols) {
  well_names <- c()
  for (row in LETTERS[1:rows]) {
    for (col in 1:cols) {
      well_names <- c(well_names, paste0(row, col))
    }
  }
  return(well_names)
}

# Set the dimensions of the plot grid
grid_rows <- 8
grid_cols <- 12

# Generate well names based on the grid dimensions
well_names <- generate_well_names(grid_rows, grid_cols)

# Create an empty dataframe with two columns: "Well" and "Condition"
design_matrix <- reactiveValues(
  data = data.frame(
    Well = well_names,
    Condition = rep(NA, length(well_names)),
    stringsAsFactors = FALSE
  )
)

# Add well names to the "Well" column in the design_matrix
design_matrix$Well <- well_names

# Define a reactiveVal to store OD data
od_matrix <- reactiveVal(NULL)

# Define a reactiveVal to store Fluorescence data
flu_matrix <- reactiveVal(NULL)

# Create a reactiveVal to store the blank_matrix dataframe
blank_matrix <- reactiveVal(NULL)

# Define a reactiveVal to store clicked points and colors
rv <- reactiveValues(
  clicked_points = data.frame(x = numeric(0), y = character(0), y_inv = character(0), id = character(0)),
  colors = rep("black", length(well_names)),
  showLegend = FALSE
)

# Define the UI for the Shiny app
ui <- fluidPage(
  titlePanel("Formatting raw plate reader data"),
  shinyjs::useShinyjs(),  # Enable the use of shinyjs
  sidebarLayout(
    sidebarPanel(
      selectInput("plateReader", "Plate reader", c("TECAN M200 Pro", "EnSight", "SpectraMax M2e", "TECAN MPlex")),
      selectInput("dataFormat", "Data format", c("Block-shape", "Wide-shape", "Tidy-shape")),
      fileInput("odFile", "OD", accept = c(".csv", ".txt", ".xlsx", ".xls")),
      fileInput("fluorescenceFile", "Fluorescence", accept = c(".csv", ".txt", ".xlsx", ".xls")),
      textInput("experimentalCondition", "Experimental Condition", value = "0"),
      actionButton("restart", "Restart"),
      actionButton("saveDesign", "Save design"),
      br(),  # Added blank space
      br(),  # Added blank space
      br(),  # Added blank space
      textInput("fileName", "File Name", value = "biosensor")  # Added text input for file name
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Design", 
                 h3("Design Matrix"),  # Title added here
                 plotOutput("pointPlot", click = "plot_click"),  # Updated variable name
                 conditionalPanel(
                   condition = "input.showLegend",
                   uiOutput("legendToggle")
                 )
        ),
        tabPanel("Data", 
                 h3("Raw plate reader data"),
                 verbatimTextOutput("odData"),
                 verbatimTextOutput("fluData")
                 # Add content for the second tab as needed
                 # ...
        ),
        tabPanel("Blank", 
                 h3("Blank matrix"),
                 fluidRow(
                   column(6,  # Adjust column width as needed
                          selectInput("assignBlanks", "Assigning blanks", choices = NULL)),
                   column(6,  # Adjust column width as needed
                          actionButton("exportBlanks", "Export")
                   )
                 ),
                 plotOutput("blank_grid", click = "blank_grid_click"),
                 verbatimTextOutput("blank_matrix_output")  # Display the blank_matrix dataframe
        )
      )
    )
  )
)

# Define the server logic for the Shiny app
server <- function(input, output, session) {
  
  # Create a reactive value to track whether "Save design" button is clicked
  saveDesignClicked <- reactiveVal(FALSE)
  
  # Create a data frame with point coordinates
  points_data <- expand.grid(
    x = seq(1, grid_cols),
    y = LETTERS[1:grid_rows]
  )
  
  # Plot the grid of points
  output$pointPlot <- renderPlot({
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, grid_cols + 1), ylim = c(0, grid_rows + 1), xaxt = "n", yaxt = "n")
    axis(1, at = 1:grid_cols, labels = 1:grid_cols)
    axis(2, at = 1:grid_rows, labels = LETTERS[grid_rows:1], las = 2)
    abline(h = 1:grid_rows, col = "lightgray", lty = 2)
    abline(v = 1:grid_cols, col = "lightgray", lty = 2)
    
    # Update colors based on unique ids
    unique_ids <- unique(rv$clicked_points$id)
    colors <- rainbow(length(unique_ids))
    
    for (i in seq_along(unique_ids)) {
      subset_points <- rv$clicked_points[rv$clicked_points$id == unique_ids[i], ]
      points_to_color <- paste(subset_points$x, subset_points$y, sep = "_")
      rv$colors[paste(points_data$x, points_data$y, sep = "_") %in% points_to_color] <- colors[i]
    }
    
    points(points_data$x, as.numeric(factor(points_data$y)), pch = 19, cex = 2, col = rv$colors)
    # Add legend if showLegend is TRUE
    if (rv$showLegend) {
      legend("top", legend = unique_ids, fill = colors, title = "Experimental Conditions", horiz = TRUE)
    }
  })
  
  
  # Plot for the blank grid in the Blank tab
  output$blank_grid <- renderPlot({
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, grid_cols + 1), ylim = c(0, grid_rows + 1), xaxt = "n", yaxt = "n")
    axis(1, at = 1:grid_cols, labels = 1:grid_cols)
    axis(2, at = 1:grid_rows, labels = LETTERS[grid_rows:1], las = 2)
    abline(h = 1:grid_rows, col = "lightgray", lty = 2)
    abline(v = 1:grid_cols, col = "lightgray", lty = 2)
    
    # Identify coordinates with id equal to 0 or "blank"
    blank_coordinates <- rv$clicked_points[rv$clicked_points$id %in% c(0, "blank"), ]
    
    # Initialize colors vector to black
    colors <- rep("black", nrow(points_data))
    
    # Change color of points to grey for identified coordinates
    if (!is.null(blank_coordinates) && nrow(blank_coordinates) > 0) {
      points_to_color <- paste(blank_coordinates$x, blank_coordinates$y, sep = "_")
      
      # Update colors for matching points to grey
      matching_indices <- which(paste(points_data$x, points_data$y, sep = "_") %in% points_to_color)
      colors[matching_indices] <- "grey"
    }
    
    # Change color of points
    points(points_data$x, as.numeric(factor(points_data$y)), pch = 19, cex = 2, col = colors)
  })
  
  
  
  
  # Render the blank_matrix dataframe in the "Blank" tab
  output$blank_matrix_output <- renderPrint({
    if (!is.null(blank_matrix())) {
      return(blank_matrix())
    } else {
      return("Blank matrix not available")
    }
  })
  
  # Update clicked points and colors when a point is clicked
  observeEvent(input$plot_click, {
    click_x <- round(input$plot_click$x)
    click_y <- LETTERS[round(input$plot_click$y)]
    
    existing_points <- rv$clicked_points[rv$clicked_points$x == click_x & rv$clicked_points$y == click_y, ]
    
    if (nrow(existing_points) > 0) {
      rv$clicked_points <- rv$clicked_points[!(rv$clicked_points$x == click_x & rv$clicked_points$y == click_y), ]
      rv$colors[paste(points_data$x, points_data$y, sep = "_") == paste(click_x, click_y, sep = "_")] <- "black"
    } else {
      condition <- input$experimentalCondition
      id <- condition
      inverted_complement_y <- LETTERS[8 - (match(click_y, LETTERS) - 1)]
      new_point <- data.frame(x = click_x, y = click_y, y_inv = inverted_complement_y, id = id)
      rv$clicked_points <- rbind(rv$clicked_points, new_point)
      colors <- rainbow(nrow(rv$clicked_points))
      rv$colors[paste(points_data$x, points_data$y, sep = "_") == paste(click_x, click_y, sep = "_")] <- colors[nrow(rv$clicked_points)]
    }
  })
  
  # Restart button
  observeEvent(input$restart, {
    rv$clicked_points <- data.frame(x = numeric(0), y = character(0), y_inv = character(0), id = character(0))
    rv$colors <- rep("black", length(well_names))
    rv$showLegend <- FALSE
  })
  
  # Dynamically render the legend toggle button
  output$legendToggle <- renderUI({
    actionButton("legendToggle", if (rv$showLegend) "Hide Legend" else "Show Legend")
  })
  
  # Save design button
  observeEvent(input$saveDesign, {
    rv$colors <- rep("white", length(well_names))
    unique_ids <- unique(rv$clicked_points$id)
    colors <- rainbow(length(unique_ids))
    
    # Create blank_matrix dataframe
    blank_matrix_data <- data.frame(
      Wells = rep(NA, length(setdiff(unique_ids, c(0, "blank")))),
      Condition = setdiff(unique_ids, c(0, "blank"))
    )
    
    rv$showLegend <- TRUE
    
    for (i in seq(nrow(rv$clicked_points))) {
      well_name <- paste(rv$clicked_points$y_inv[i], rv$clicked_points$x[i], sep = "")
      row_index <- which(design_matrix$data$Well == well_name)
      
      if (length(row_index) > 0) {
        design_matrix$data$Condition[row_index] <- rv$clicked_points$id[i]
      
      }
    }
    
    # Update blank_matrix
    blank_matrix(blank_matrix_data)
  })
  
  
  
  # Define a function to handle y-coordinate inversion
invert_y_coordinate <- function(y, total_rows) {
  return(toupper(letters[total_rows + 1 - y]))
}

# Observer for the "blank grid" click
observeEvent(input$blank_grid_click, {
  click_x <- round(input$blank_grid_click$x)
  click_y <- invert_y_coordinate(round(input$blank_grid_click$y), grid_rows)  # Adjust for the inverted y-axis
  
  # Find the clicked well in the blank_matrix_data
  clicked_well <- paste(click_y, click_x, sep = "")
  
  blank_matrix_data <- blank_matrix()
  
  if (!is.null(blank_matrix_data) && nrow(blank_matrix_data) > 0) {
    # Check if the current value of the "Selection box" is in the "Condition" column
    condition_value <- input$assignBlanks
    row_index <- which(blank_matrix_data$Condition == condition_value)
    
    if (length(row_index) > 0) {
      # Update the corresponding "Well" value with the clicked coordinate
      blank_matrix_data$Wells[row_index] <- clicked_well
      
      # Update the blank_matrix
      blank_matrix(blank_matrix_data)
    }
  }
})

  # Inside the observeEvent for "odFile" input
  observeEvent(input$odFile, {
    if (!is.null(input$odFile)) {
      od_matrix(block_shape_processing(input$odFile$datapath))
    } else {
      od_matrix(NULL)
    }
  })
  
  # Inside the observeEvent for "fluorescenceFile" input
  observeEvent(input$fluorescenceFile, {
    if (!is.null(input$fluorescenceFile)) {
      flu_matrix(block_shape_processing(input$fluorescenceFile$datapath))
    } else {
      flu_matrix(NULL)
    }
  })
  
  # Update the UI to display the data frames or the message
  output$odData <- renderPrint({
    if (!is.null(od_matrix())) {
      return(od_matrix())
    } else {
      return("File non existent")
    }
  })
  
  output$fluData <- renderPrint({
    if (!is.null(flu_matrix())) {
      return(flu_matrix())
    } else {
      return("File non existent")
    }
  })
  
  # Update choices for "Assigning blanks" selectInput
  observe({
    unique_ids <- unique(rv$clicked_points$id)
    valid_ids <- unique_ids[!(unique_ids %in% c(0, "blank"))]
    updateSelectInput(session, "assignBlanks", choices = valid_ids)
  })
  
  # Observer for "Export" button
  observeEvent(input$exportBlanks, {
    # Save the data frames to separate CSV files
    write.csv(od_matrix(), file = "od_matrix.csv", row.names = FALSE)
    write.csv(flu_matrix(), file = "flu_matrix.csv", row.names = FALSE)
    write.csv(design_matrix$data, file = "design_matrix.csv", row.names = FALSE)
    
    blank_matrix_data <- blank_matrix()
    if (!is.null(blank_matrix_data)) {
      write.csv(blank_matrix_data, file = "blank_matrix.csv", row.names = FALSE)
    } else {
      print("No data available for blank_matrix.")
    }
  })
  
}

# Run the application
shinyApp(ui, server)