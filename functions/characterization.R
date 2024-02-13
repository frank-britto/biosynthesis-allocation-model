
titration <- function(name,
                      time_limits) {
  
  # ########################################################################################
  # titration function:
  # 
  # Description: scatter plot of phi values with respect to the condition
  #               
  #
  # Inputs: 
  #         name:         data-frame obtained from 'cleaning_data.R' <type:tibble>
  #         time_limits:  time interval from which the hill function will be constructed
  #                       <type: list>
  #
  # Output: 
  #         max_mu:       matrix with the phi values or slopes with their corresponding
  #                       condition value, according to the chosen method <type: tibble>
  #
  # #########################################################################################
  
  # Obtain the leaky expression level (baseline)
  baseline <- name %>%
    filter(time > time_limits[1]) %>%  
    filter(time < time_limits[2]) %>%
    filter(Condition == min(Condition))
  
  # Obtain data within time limits
  tab <- name %>%
    filter(time > time_limits[1]) %>%  
    filter(time < time_limits[2])
  
  # Substract the baseline expression
  tab$phi <- tab$phi - baseline$phi
  
  tab <- tab %>%
    filter(Condition > min(Condition))
  
  # Obtain the maximum value of heterologous fraction
  max_mu <- tab %>%
    group_by(Condition) %>%
    filter(phi == max(phi)) %>%
    select(phi)
  
  return(max_mu)
  
}

dose_response <- function(max_mu, 
                          initial) {
  
  # ########################################################################################
  # dose_response:
  # 
  # Description: non-linear model for fitting data
  #
  # Inputs: 
  #         max_mu:       data-frame obtained from 'titration.R' <type:tibble>
  #         max_ind:      maximun concentration of Condition <unit: µM>
  #         initial:      initial guesses for non-linear fitting. <type: list>
  #                       For example, initial <- list(A=8000, ki=50, n=2)
  #
  # Output: 
  #         coeff:        coefficients from the linear regression <type: vector>
  #         model_fit:    non-linear regression based on the hill function <type: list>
  #
  # #######################################################################################
  
  # Estimate the hill function coefficients
  model <- max_mu %>%
    nls(phi~H*(Condition^n / (Condition^n + ki^n)), 
        data = ., start = initial,
        algorithm="port", lower=c(0,0,1))
  
  parameters <- model %>% 
    summary() %>% .$parameters %>% 
    as_tibble() %>%
    mutate(parameter =c("H", "ki", "n"), .before = Estimate) %>%
    pull(Estimate)
  
  # Saving hill function coefficients
  H <- coef(model)["H"]
  ki <- coef(model)["ki"]
  n <- coef(model)["n"]
  
  coeff <- coef(model)
  
  # Calculate fitting errors
  iptg <- 0:max(max_mu$Condition)
  fi_hat <- predict(model, newdata = tibble(Condition = iptg))
  
  parameter_standard_errors <- sqrt(diag(vcov(model)))
  
  prediction_error <- sqrt(
    ((iptg)^n / ((iptg)^n + ki^n))^2 * parameter_standard_errors["H"]^2 +
      (H * n * ki^(n-1) * (iptg)^n / ((iptg)^n + ki^n)^2 )^2 * parameter_standard_errors["ki"]^2 +
      (H * ki^n * log(iptg / ki) * (iptg)^n / ((iptg)^n + ki^n)^2 )^2 * parameter_standard_errors["n"]^2
  )
  
  model_fit <- tibble(Condition = 0:max(max_mu$Condition), fi = fi_hat,
                      # Calculate lower and upper bounds for the confidence intervals
                      lower = fi_hat - qnorm(0.975) * prediction_error, 
                      upper = fi_hat + qnorm(0.975) * prediction_error, 
                      error = prediction_error)
  
  return(list(model_fit, model))
  
}


dose_response_visualization <- function(name, time_limits, initial_guess, flag = FALSE) {
  
  # ########################################################################################
  # dose_response:
  # 
  # Description: non-linear model for fitting data
  #
  # Inputs: inputs are the same as those used for the titration and dose_response functions
  #         switch flag = TRUE if the user wants to plot the fitted curve
  #
  # Output: 
  #         dose response curve <type: ggplot2>
  #
  # #######################################################################################
  
  experiment <- titration(name = name,
                          time_limits = time_limits
                          )
  
  dose_response_curve <- experiment %>%
                          ggplot(data = ., aes(Condition, phi))+
                          geom_point()+
                          xlab('Condition (µM)')+
                          ylab(expression(FU~OD^-1))+
                          ggtitle("Hill function")+
                          theme(plot.title = element_text(size = 20, face = "bold"))
    
    
    
  if (flag == TRUE) {
    
    output <- dose_response(experiment,
                            initial = initial_guess)
    
    dose_response_curve <- dose_response_curve +
      geom_line(data = output[[1]], aes(Condition, fi))
    
    return(list(dose_response_curve, output[[2]]))
    
  }
  
  else {
    
    return(list(dose_response_curve, NULL))
    
  }
  
}

slope_based_visualization <- function(linear_df, initial_guess, flag = TRUE) {
  
  # Change column named "Slope" for "phi"
  colnames(linear_df)[colnames(linear_df) == "Slope"] <- "phi"
  
  dose_response_curve <- linear_df %>%
                            ggplot(data = ., aes(Condition, phi))+
                            geom_point()+
                            xlab('Condition (µM)')+
                            ylab(expression(FU~OD^-1))+
                            ggtitle("Hill function")+
                            theme(plot.title = element_text(size = 20, face = "bold"))
  
                          
  if (flag == TRUE) {
    
    output <- dose_response(linear_df,
                            initial = initial_guess)
    
    dose_response_curve <- dose_response_curve +
      geom_line(data = output[[1]], aes(Condition, fi))
    
    return(list(dose_response_curve, output[[2]]))
    
  }
  
  else {
    
    return(list(dose_response_curve, NULL))
    
  }
  
  
}

annotations_hook_plot <- function(name,
                                  end) {
  
  # ######################################################################################
  # annotations_hook_plot function: 
  # 
  # Description: locations the time points where where the heterologous fraction reaches 
  # equilibrium (maximum) and the maximum growth rate
  #
  # Inputs: 
  #         name:         data-frame obtained from 'cleaning_data.R' <type:tibble>
  #         time_limits:  lower and upper time bounds <type:float>
  #                 
  # Output: 
  #         tab_Hwa:      locations of maximum growth rate
  #         tab_Vacc:     locations of maximum heterologous fraction
  #
  # #######################################################################################
  
  # Obtaining the baseline expression
  bl <- name  %>%
    filter(Condition == min(Condition)) %>% 
    select(time, phi) %>% 
    rename(bl = phi) %>% 
    ungroup()
  
  # Substracting the baseline
  tidy <- name %>% 
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
              gr_index=which.max(growth_rate))
  
  aux <- name %>% filter(Condition > min(Condition))
  a <- unique(aux$Condition)
  b <- indexes$gr_index
  c <- indexes$phi_index
  
  # Dataframe containing the location of maximum growth rate
  tab_Hwa <- sapply(2:length(a),function(i){
    tab <- tidy %>% filter(Condition==a[i]) %>%
      select(production_rate, growth_rate) %>% .[b[i],]
    tab <- bind_cols(tab, a[i])
  }) %>% t() %>% as_tibble %>% 
    unnest(cols = c(`...3`, production_rate, growth_rate)) %>% 
    rename(Condition = `...3`) %>% relocate(Condition)
  
  # Dataframe containing phi at equilibrium
  tab_Vacc<- sapply(2:length(a),function(i){
    tab <- tidy %>% filter(Condition==a[i]) %>%
      select(production_rate, growth_rate) %>% .[c[i],]
    tab <- bind_cols(tab, a[i])
  }) %>% t() %>% as_tibble %>% 
    unnest(cols = c(`...3`, production_rate, growth_rate)) %>% 
    rename(Condition = `...3`) %>% relocate(Condition) 
  
  return(list(tab_Hwa, tab_Vacc))
  
}

hook_plot <- function(name, time_limits) {
  
  # ######################################################################################
  # hook_plot function: 
  # 
  # Description: generate the production rate vs. growth rate plot, and locate the
  #              time where the heterologous fraction reaches equilibrium
  #
  # Inputs: 
  #         name:         data-frame obtained from 'cleaning_data.R' <type:tibble>
  #         time_limits:  lower and upper time bounds <type:float>
  #                 
  # Output: 
  #         ggplot2 object with the annotated hook plot
  #
  # #######################################################################################
  
  locations <- annotations_hook_plot(name, time_limits[[2]])
  
  # Generating the hook plot
  name %>% 
    filter(time > time_limits[1] & time < time_limits[2]) %>%
    mutate(Condition = as.factor(Condition)) %>%
    ggplot(data = ., aes(growth_rate, production_rate/1000, color = Condition)) +
    geom_path(size = 1.5)+
    theme_classic()+
    scale_colour_grey()+
    geom_point(data = locations[[2]], aes(growth_rate, production_rate/1000), 
               color = 'green', size = 5, shape=22, stroke = 2)+
    geom_point(data = locations[[1]], aes(growth_rate, production_rate/1000), 
               color = 'red3', size = 5, shape = 21, stroke = 2)+
    xlab(expression(h^-1))+
    ylab(expression(10^3~FU~OD^-1~h^-1))+
    ggtitle("Hook plot")+
    theme(plot.title = element_text(size = 15, face = "bold"))
}

linear_detection <- function(name, window_size = 20) {
  
  calculate_r_squared <- function(x, y) {
    lm_model <- lm(y ~ x)
    summary(lm_model)$r.squared
  }
  
  calculate_highest_r_squared_interval_sliding_window <- function(production_data, growth_data, window_size) {
    n_points <- nrow(growth_data)
    max_r_squared <- -Inf
    max_interval <- c(0, 0)
    max_slope <- NA_real_
    
    for (i in 1:(n_points - window_size + 1)) {
      # Define the current window
      current_window <- growth_data[i:(i + window_size - 1), ]
      
      # Filter production_rate based on the current window
      filtered_production_data <- production_data %>%
        filter(time >= min(current_window$time), time <= max(current_window$time))
      
      # Perform linear regression
      lm_model <- lm(filtered_production_data$values ~ current_window$values)
      
      # Check if the slope is positive
      if (lm_model$coefficients[2] > 0) {
        # Calculate R-squared
        current_r_squared <- summary(lm_model)$r.squared
        
        # Update max R-squared and corresponding interval
        if (current_r_squared > max_r_squared) {
          max_r_squared <- current_r_squared
          max_interval <- c(min(current_window$time), max(current_window$time))
          max_slope <- lm_model$coefficients[2]
        }
      }
    }
    
    return(list(interval = max_interval, slope = max_slope))
  }
  
  gr <- tibble(
    "time" = name$time,
    "values" = name$growth_rate
  )
  
  pr <- tibble(
    "time" = name$time,
    "values" = name$production_rate
  )
  
  return(calculate_highest_r_squared_interval_sliding_window(pr, gr, window_size))
  
}

generate_linear_detection_results <- function(input_df) {
  # Extract unique values from the "Condition" column
  conditions <- unique(input_df$Condition)
  
  # Create a new dataframe to store the results
  result_df <- tibble(
    "Condition" = conditions,
    "Minimum" = rep(NA_real_, length(conditions)),
    "Maximum" = rep(NA_real_, length(conditions)),
    "Slope" = rep(NA_real_, length(conditions))
  )
  
  # Iterate over each condition
  for (i in seq_along(conditions)) {
    # Subset the dataframe for the current condition
    subset_df <- input_df %>% filter(Condition == conditions[[i]])
    
    # Apply the linear_detection function
    linear_detection_result <- linear_detection(subset_df)
    
    # Directly assign the values to the specific row indices
    result_df$Minimum[i] <- linear_detection_result$interval[1]
    result_df$Maximum[i] <- linear_detection_result$interval[2]
    result_df$Slope[i] <- linear_detection_result$slope
  }
  
  return(result_df)
}


visualize_manual_selection <- function(name, t1, t2, condition) {
  # Fit a linear model
  lm_model <- lm(production_rate/1000 ~ growth_rate, data = name %>% 
                   filter(Condition == condition, time > t1, time < t2))
  
  # Access the R-squared value
  r_squared <- summary(lm_model)$r.squared
  
  # Plot with time filter
  plot_with_time_filter <- name %>%
    filter(Condition == condition, time > t1, time < t2) %>%
    ggplot(aes(growth_rate, production_rate/1000)) +
    geom_point(size = 3) +
    geom_abline(slope = coef(lm_model)[2], intercept = coef(lm_model)[1], color = "red", size = 1.5) +  # Bigger linewidth
    ggtitle(paste("R^2 =", sprintf("%.3f", r_squared))) +
    theme_classic() +
    xlab(expression(h^-1)) +
    ylab(expression(10^3~FU~OD^-1~h^-1))
  
  
  
  # Plot without time filter
  plot_without_time_filter <- name %>%
    filter(Condition == condition) %>%
    ggplot(aes(growth_rate, production_rate/1000)) +
    geom_path(size = 3) +
    geom_abline(slope = coef(lm_model)[2], intercept = coef(lm_model)[1], color = "red", size = 1.5) +  # Bigger linewidth
    theme_classic() +
    xlab(expression(h^-1)) +
    ylab(expression(10^3~FU~OD^-1~h^-1))
  
  return(list(plot_with_time_filter, plot_without_time_filter))
}
