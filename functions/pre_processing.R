# #########################################################################
# CUSTOM FUNCTIONS FOR NON-FORMATTED DATASETS
# #########################################################################

excel2tibble <- function(name) {
  
  # ######################################################################## 
  # Description: process data from TECAN* from a SINGLE strain. It           
  # assumes that we are working inside the folder where the data is located  
  #                                                                         
  # Inputs:                                                                  
  #         name:         name of the Excel file's path <type:string>         
  #                                                                          
  # Output:                                                                  
  #         tidy:         output tidy dataframe <type:tibble>                       
  #                                                                          
  # Comments: this is the type of data used for making all plots.            
  # ######################################################################## 
  
  # Vector to store the different inductor concentrations
  iptg <- c()
  
  # Obtain the strain name from the file
  file_name <- basename(name)  
  components <- unlist(strsplit(file_name, "_"))  
  strain_name <- components[1]
  
  # Processing each Excel tab
  sheets <- strsplit(excel_sheets(name), split = " ")
  for (j in seq_along(1:length(sheets))){
    iptg <- append(iptg, as.numeric(sheets[[j]][1])*1000)
  }
  
  for (ind in seq_along(1:length(iptg))) {
    
    # Obtaining the index of the empty column
    idx <- NA
    tab <- read_excel(name, sheet = excel_sheets(name)[ind])
    
    for (i in seq_along(colnames(tab))){
      if (substring(colnames(tab)[i],1,1) == '.') {
        idx <- i
        break
      }
    }
    
    # Get the OD and fluorescence data
    t <- tab %>% select(colnames(tab)[1])
    od <- tab %>% select(colnames(tab)[2:as.integer(idx-1)]) %>% mutate_all(~.-0.1)
    flu <- tab %>% select(colnames(tab)[as.integer(idx+2): length(colnames(tab))])
    
    # Calculate the derived variables
    phi <- flu / od
    diff_od <- apply(od, 2, diff)
    diff_flu <- apply(flu, 2, diff)
    diff_t <- apply(t, 2, diff)
    gr <- diff_od/(od*diff_t)
    pr <- diff_flu/(od*diff_t)
    
    # Applying mean
    od <- od %>% apply(1,mean)
    flu <- flu %>% apply(1,mean)
    phi <- phi %>% apply(1,mean)
    gr <- gr %>% apply(1,mean)
    pr <- pr %>% apply(1,mean)
    
    # Generate matrix with all this data
    base <- bind_cols(t, od, flu, phi, gr, pr) %>% 
      mutate(inductor = iptg[ind]) %>%
      set_names(c("time", "od", "flu", "phi", "gr", "pr", "inductor"))
    
    if (ind == 1){
      bases <- base
    } 
    else{
      bases <- bind_rows(bases, base)
    }
    
  }

  tidy <- tibble(
    time = bases$time,
    od = bases$od,
    flu = bases$flu,
    phi = bases$phi,
    growth_rate = bases$gr,
    production_rate = bases$pr,
    Condition = bases$inductor,
    strain = strain_name
  )
  
  return(tidy)
  
}

# #########################################################################
# FUNCTIONS FOR FORMATTED DATASETS
# #########################################################################

block_shape_processing <- function(name) {
  
  # ######################################################################## 
  # Description:  block-shape to wide-shape format  
  #                                                                         
  # Inputs:                                                                  
  #               name:         name of the Excel file's path <type:string>         
  #                                                                          
  # Output:                                                                  
  #               raw_data:     output tidy dataframe <type:tibble>                                  
  # ########################################################################
  
  # Read the Excel file in the block-shape format
  tab <- read_excel(name, col_names = FALSE)
  
  # Generate the time interval vector (TECAN format: each number is followed by an "s")
  time <- tab$...1[grep("s$", tab$...1)]
  time <- as.numeric(sub("s$", "", time, ignore.case = TRUE))
  time <- time/3600 # convert to hours
  
  # Generate the column names for the wide-shape dataframe
  letters_vector <- LETTERS[1:8]
  numbers_vector <- 1:12
  
  # Create a vector to store the combinations
  well_columns <- c()
  
  # Loop over letters and numbers to generate combinations
  for (letter in letters_vector) {
    for (number in numbers_vector) {
      # Concatenate the combinations and append to the vector
      well_columns <- c(well_columns, paste0(letter, number))
    }
  }
  
  # Generating the wide-shape dataframe "scaffold"
  raw_data <- data.frame(matrix(NA, nrow = length(time), ncol = length(well_columns)))
  
  # Set column names to well_columns
  colnames(raw_data) <- well_columns
  
  # Dropping the first column
  tab <- select(tab, -...1)
  
  # Iterating over every "snapshot"
  snapshot <- 12 - 1 
  num_rows <- nrow(tab) 
  num_iterations <- num_rows %/% snapshot
  
  for (i in seq(1, num_iterations, by = 1)) {
    
    start_row <- (i - 1) * snapshot + 1
    end_row <- min(i * snapshot, num_rows)
    
    chunk <- tab[start_row:end_row, ]

    # Drop the first three rows
    chunk <- chunk[-c(1, 2, 3), ]
    
    # Replace the ith row of raw_data with reshaped vector
    raw_data[i, ] <- as.vector(t(as.matrix(chunk)))
  }
  
  # Adding the "time" interval
  raw_data$time <- time
  raw_data <- raw_data[, c(ncol(raw_data), 1:(ncol(raw_data)-1))]
  
  return(raw_data)
  
}

raw2tidy <- function(od_matrix, flu_matrix, design_matrix, blank_matrix) {
  
  # Not finished
  
  # Find matching Wells values with the same Condition
  matching_wells <- design_matrix %>%
    group_by(Condition) %>%
    pull(Well)
  
  # Extract corresponding columns from od_matrix
  selected_columns <- od_matrix[, matching_wells]
  
  # Calculate row-wise mean
  mean_values <- rowMeans(selected_columns)
  
  # Create a new dataframe
  result_dataframe <- data.frame(Mean_Values = mean_values, Condition = design_matrix$Condition)
  
  return(result_dataframe)
  
}
