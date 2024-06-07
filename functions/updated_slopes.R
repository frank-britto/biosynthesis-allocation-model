
smooth_filter <- function(df, window_size = 5) {
  
  # ######################################################################################
  # smooth_filter function: 
  # 
  # Description: rolling window algorithm for smoothing the "hook" plot
  #
  # Inputs: 
  #         df:           data-frame obtained from 'raw2tidy' <type:tibble>
  #         window_size:  size of the window <type:integer>
  #                 
  # Output: 
  #         df:           same data-frame with new columns:
  #                       "gr": re-computing of growth rate with smoothing
  #                       "pr": re-computing of production rate with smoothing
  #
  # #######################################################################################
  
  # Initialize a vector to store the growth rates
  growth_rates <- numeric(nrow(df))
  production_rates <- numeric(nrow(df))
  
  # Perform rolling window calculation
  for (i in 1:(nrow(df) - window_size + 1)) {
    
    # Subset the data within the window
    window_data <- df[i:(i + window_size - 1), ]
    
    # Fit linear regression model with "od" as dependent variable
    model1 <- lm(od ~ time, data = window_data)
    model2 <- lm(flu ~ time, data = window_data)
    
    # Extract slope coefficient (growth rate)
    growth_rate <- coef(model1)[2]  
    production_rate <- coef(model2)[2]
    
    # Assign the growth rate to the corresponding index
    growth_rates[i] <- growth_rate
    production_rates[i] <- production_rate
  }
  
  # Pad with NA values for incomplete windows at the end
  growth_rates <- c(growth_rates, rep(NA, nrow(df) - length(growth_rates)))
  production_rates <- c(production_rates, rep(NA, nrow(df) - length(production_rates)))
  
  # Add growth_rates to the dataframe
  df$gr <- growth_rates
  df$pr <- production_rates
  
  return(df)
  
}

rates_plot <- function(data,
                      linear_section_matrix) {
  
  hplot <- data %>%
    mutate(Condition = as.factor(Condition)) %>%
    ggplot(aes(gr, pr, color = time)) +  
    geom_path(size = 1.5) +
    xlab(expression(μ~(h^-1))) +
    ylab(expression(Φ[H]~(10^3~FU)~(OD^-1)~(h^-1)))+
    scale_color_gradient(low = "lightgrey", high = "black") +
    facet_wrap(~Condition) 
  
  linear_regression <- hplot +
    geom_smooth(data = . %>% filter(time >= linear_section_matrix$Minimum[1] & time <= linear_section_matrix$Maximum[1]),
                method = "lm", se = FALSE, color = "#1da9e0")
  
  return(linear_regression)
  
}

annotations_hook_plot <- function(data,
                                  end) {
  
  # ######################################################################################
  # annotations_hook_plot function: 
  # 
  # Description: locations the time points where where the heterologous fraction reaches 
  # equilibrium (maximum) and the maximum growth rate
  #
  # Inputs: 
  #         data:         data-frame obtained from 'cleaning_data.R' <type:tibble>
  #         time_limits:  lower and upper time bounds <type:float>
  #                 
  # Output: 
  #         tab_Hwa:      locations of maximum growth rate
  #         tab_Vacc:     locations of maximum heterologous fraction
  #
  # #######################################################################################
  
  # Obtaining the baseline expression from the whole dataset
  bl <- data  %>%
    filter(Condition == min(Condition)) %>% 
    select(time, phi) %>% 
    rename(bl = phi) %>% 
    ungroup()
  
  # Subtracting the baseline
  tidy <- data %>% 
    left_join(bl, by = c("time")) %>%
    filter(Condition > min(Condition)) %>% 
    mutate(phi = phi - bl)%>%
    mutate(Condition = as.numeric(Condition))%>%
    arrange(Condition)%>%
    mutate(Condition = as.factor(Condition)) 
  
  # Obtaining the indexes for maximum growth rate and phi at equilibrium
  indexes <- tidy %>% filter(between(time, 0, end)) %>% 
    group_by(Condition) %>%
    summarize(phi_index=which.max(phi),
              gr_index=which.max(gr)) # growth_rate
  
  # Creating auxiliary variables
  aux <- data %>% filter(Condition > min(Condition))
  a <- unique(aux$Condition)
  b <- indexes$gr_index
  c <- indexes$phi_index
  
  # Dataframe containing the location of maximum growth rate
  tab_Hwa <- sapply(2:length(a),function(i){
    tab <- tidy %>% filter(Condition==a[i]) %>%
      select(pr, gr) %>% .[b[i],]
    tab <- bind_cols(tab, a[i])
  }) %>% t() %>% as_tibble %>% 
    unnest(cols = c(`...3`, pr, gr)) %>% 
    rename(Condition = `...3`) %>% relocate(Condition)
  
  # Dataframe containing phi at equilibrium
  tab_Vacc<- sapply(2:length(a),function(i){
    tab <- tidy %>% filter(Condition==a[i]) %>%
      select(pr, gr) %>% .[c[i],]
    tab <- bind_cols(tab, a[i])
  }) %>% t() %>% as_tibble %>% 
    unnest(cols = c(`...3`, pr, gr)) %>% 
    rename(Condition = `...3`) %>% relocate(Condition) 
  
  return(list(tab_Hwa, tab_Vacc))
  
}

hook_plot <- function(data, time_limits, highlight = FALSE, condition = NaN) {
  
  # ######################################################################################
  # hook_plot function: 
  # 
  # Description: generate the production rate vs. growth rate plot, and locate the
  #              time where the heterologous fraction reaches equilibrium
  #
  # Inputs: 
  #         data:         data-frame obtained from 'cleaning_data.R' <type:tibble>
  #         time_limits:  lower and upper time bounds <type:float>
  #                 
  # Output: 
  #         ggplot2 object with the annotated hook plot
  #
  # #######################################################################################
  
  locations <- annotations_hook_plot(data, time_limits[[2]])
  
  hplot <- data %>%
    mutate(Condition = as.factor(Condition)) %>%
    ggplot(aes(gr, pr, color = Condition)) +  
    xlab(expression(μ~(h^-1))) +
    ylab(expression(Φ[H]~(10^3~FU)~(OD^-1)~(h^-1)))
  
  # Highlight one plot
  if (highlight == TRUE) {

    hplot <- hplot +
      geom_path(size = 1.5, alpha = ifelse(data$Condition != condition, 0.2, 1)) +
      geom_point(data = locations[[2]] %>% filter(Condition == condition), aes(gr, pr), 
                 color = 'black', size = 5, shape=22, stroke = 2)+
      geom_point(data = locations[[1]] %>% filter(Condition == condition), aes(gr, pr), 
                 color = 'red3', size = 5, shape = 21, stroke = 2)
  }
  
  else {
    
    # Generating the hook plot
    hplot <- hplot +
      geom_path(size = 1.5)+
      geom_point(data = locations[[2]], aes(gr, pr), 
                 color = 'black', size = 5, shape=22, stroke = 2)+
      geom_point(data = locations[[1]], aes(gr, pr), 
                 color = 'red3', size = 5, shape = 21, stroke = 2)
    
  }
  
  return(hplot)
  
  
}

preliminary_slopes <- function(data, time_interval) { 
  
  # ######################################################################################
  # preliminary_slopes function: 
  # 
  # Description: generate a dataframe with the linear interval with the highest R2
  #              value calculated by the linear_detection function
  #
  # Inputs: 
  #         data:       data-frame obtained from 'cleaning_data.R' <type:tibble>
  #         time_interval:  time interval selected in the "Data summary" section
  #                 
  # Output: 
  #         dataframe with the linear interval with the highest R2
  #
  # #######################################################################################
  
  # Filter for the minimum condition
  data <- data %>%
    filter(Condition > min(Condition))
  
  # Extract unique values
  conditions <- unique(data$Condition)
  
  # Create a new dataframe to store the results
  result_df <- tibble(
    "Condition" = conditions,
    "Minimum" = time_interval[[1]],
    "Maximum" = time_interval[[2]],
    "Slope" = rep(NA_real_, length(conditions)),
    "R2" = rep(NA_real_, length(conditions))
  )
  
  # Iterate over each condition
  for (i in seq_along(conditions)) {
    
    # Subset the dataframe for the current condition
    subset_df <- data %>% 
      filter(Condition == conditions[[i]]) %>%
      filter(time >= time_interval[[1]] & time <= time_interval[[2]])
    
    # Perform linear regression
    lm_model <- lm(subset_df$pr ~ subset_df$gr)
    
    # Directly assign the values to the specific row indices
    result_df$Slope[i] <- lm_model$coefficients[2]
    result_df$R2[i] <- summary(lm_model)$r.squared
  }
  
  return(result_df)
}

linear_view <- function(data, condition, time_interval, t1 = NaN, t2 = NaN) {
  
  # Working with the data subset
  subset_df <- data %>% 
    filter(Condition == condition)
  
  # Defining a lower and upper bounds
  if (is.nan(t1)) {t1 <- time_interval[[1]]}
  if (is.nan(t2)) {t2 <- time_interval[[2]]}
  
  subset_data <- subset(subset_df, time >= t1 & time <= t2)
  
  # Fit a linear model
  lm_model <- lm(pr ~ gr, data = subset_data)
  
  # Access the R-squared value
  r_squared <- summary(lm_model)$r.squared
  
  # Plot with time filter
  plot_with_time_filter <- subset_data %>%
    ggplot(aes(gr, pr, color = time)) +
    geom_point(size = 3) +
    geom_abline(slope = coef(lm_model)[2], intercept = coef(lm_model)[1], color = "#1da9e0", size = 1.5) +  # Bigger linewidth
    ggtitle(bquote(R^2 == .(sprintf("%.3f", r_squared))))+
    theme_classic() + 
    scale_color_gradient(low = "lightgrey", high = "black") +
    xlab(expression(μ~(h^-1))) +
    ylab(expression(Φ[H]~(10^3~FU)~(OD^-1)~(h^-1)))
  
  # Plot without time filter
  plot_without_time_filter <- subset_df %>%
    ggplot(aes(gr, pr, color = time)) +
    geom_path(size = 3) +
    geom_abline(slope = coef(lm_model)[2], intercept = coef(lm_model)[1], color = "#1da9e0", size = 1.5) +  # Bigger linewidth
    ggtitle(paste("Condition = ", condition))+
    theme_classic() +
    scale_color_gradient(low = "lightgrey", high = "black") +
    xlab(expression(μ~(h^-1))) +
    ylab(expression(Φ[H]~(10^3~FU)~(OD^-1)~(h^-1))) #+
    #geom_point(data = subset_df, aes(x = gr, y = pr), color = "yellow", size = 3)
  
  return(list(plot_with_time_filter, plot_without_time_filter))
}

manual_adjust <- function(input_df, linear_section_matrix, condition, t1 = NaN, t2 = NaN) {
  
  # Find the row index in linear_section_matrix with the matching condition
  row_index <- which(linear_section_matrix$Condition == condition)
  
  # If the row does not exist, return the original matrix
  if (length(row_index) == 0) {
    warning("Row not found in linear_section_matrix for the given condition.")
    return(linear_section_matrix)
  }
  
  # Create a copy of the original row
  condition_row <- linear_section_matrix[row_index, ]
  
  # Update Minimum and Maximum values if t1 and/or t2 are not NaN
  if (!is.nan(t1)) {
    condition_row$Minimum <- t1
  }
  
  if (!is.nan(t2)) {
    condition_row$Maximum <- t2
  }
  
  # If either t1 or t2 is not NaN, proceed with further adjustments
  if (!is.nan(t1) || !is.nan(t2)) {
    
    # Filter input_df based on condition and time interval
    subset_df <- input_df %>% filter(time >= condition_row$Minimum & time <= condition_row$Maximum & Condition == condition)
    
    # Fit a linear model
    lm_model <- lm(pr ~ gr, data = subset_df)
    
    # Update the corresponding value of "Slope" in condition_row
    condition_row$Slope <- lm_model$coefficients[2]
    
    # Update the corresponding R-squared value
    condition_row$R2 <- summary(lm_model)$r.squared
    
    # Replace the old row in linear_section_matrix with the updated condition_row
    linear_section_matrix[row_index, ] <- condition_row
  }
  
  # Return the updated linear_section_matrix
  return(linear_section_matrix)
}

