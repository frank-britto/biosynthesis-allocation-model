
leakiness <- function(data,
                      time_interval) {
  
  # Obtain leaky expression in time
  baseline <- data %>%
    filter(time > time_interval[1]) %>%  
    filter(time < time_interval[2]) %>%
    filter(Condition == min(Condition))
  
  return(baseline)
}

titration <- function(data,
                      time_interval) {
  
  # ########################################################################################
  # titration function:
  # 
  # Description: extract max phi values for each experimental condition
  #               
  #
  # Inputs: 
  #         data:           literally our data 
  #         time_interval:  time interval defined using the "view" function (manually)
  #         baseline:       output of the "leakiness" function
  #         
  #
  # Output: 
  #         max_phi:      matrix with the phi values with their corresponding
  #                       condition value, according to the chosen method <type: tibble>
  #
  # #########################################################################################
  
  # Obtain the leaky expression level (baseline)
  baseline <- leakiness(data, time_interval)
  
  # Obtain data within time limits
  tab <- data %>%
    filter(time > time_interval[1]) %>%  
    filter(time < time_interval[2])
  
  # Substract the baseline expression
  tab$phi <- tab$phi - baseline$phi
  
  tab <- tab %>%
    filter(Condition > min(Condition))
  
  # Obtain the maximum value of heterologous fraction
  max_phi <- tab %>%
    group_by(Condition) %>%
    filter(phi == max(phi)) %>%
    select(phi)
  
  return(max_phi)
  
}

dose_response <- function(max_phi, 
                          initial_guess) {
  
  # Estimate the hill function coefficients
  model <- max_phi %>%
    nls(phi~H*(Condition^n / (Condition^n + ki^n)), 
        data = ., start = initial_guess,
        algorithm="port", lower=c(0,0,1))
  
  # Create dataframe with model parameters
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
  iptg <- 0:max(max_phi$Condition)
  fi_hat <- predict(model, newdata = tibble(Condition = iptg))
  
  parameter_standard_errors <- sqrt(diag(vcov(model)))
  
  prediction_error <- sqrt(
    ((iptg)^n / ((iptg)^n + ki^n))^2 * parameter_standard_errors["H"]^2 +
      (H * n * ki^(n-1) * (iptg)^n / ((iptg)^n + ki^n)^2 )^2 * parameter_standard_errors["ki"]^2 +
      (H * ki^n * log(iptg / ki) * (iptg)^n / ((iptg)^n + ki^n)^2 )^2 * parameter_standard_errors["n"]^2
  )
  
  model_fit <- tibble(Condition = 0:max(max_phi$Condition), fi = fi_hat,
                      # Calculate lower and upper bounds for the confidence intervals
                      lower = fi_hat - qnorm(0.975) * prediction_error, 
                      upper = fi_hat + qnorm(0.975) * prediction_error, 
                      error = prediction_error)
  
  return(list(model_fit, model))
  
}

dose_response_profiling <- function(data, time_interval, initial_guess, flag = FALSE) {
  
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
  
  # If we are working with the slope-based method
  if (any(grepl("Slope", colnames(data)))) {
    
    # Change column datad "Slope" for "phi"
    experiment <- data
    colnames(experiment)[colnames(experiment) == "Slope"] <- "phi"
  }
  
  # If we are working with the maximum method
  else {
    
    # Obtain the maximum phi values
    experiment <- titration(data = data,
                            time_interval = time_interval)
  }
  
  # Generate the dose-response curve
  dose_response_curve <- experiment %>%
    ggplot(data = ., aes(Condition, phi))+
    geom_point()+
    xlab('Condition (µM)')+
    ylab(expression(Φ[H] ~ (FU) ~ (OD^-1)))+
    ggtitle("Hill function")+
    theme(plot.title = element_text(size = 20, face = "bold"))
  
  # Fit the model
  if (flag == TRUE) {
    
    model <- dose_response(experiment,
                           initial = initial_guess)
    
    # model[[1]] = dataframe with fitting errors
    # model[[2]] = model details
    
    dose_response_curve <- dose_response_curve +
      geom_line(data = model[[1]], aes(Condition, fi))
    
    return(list(dose_response_curve, model[[2]], model[[1]]))
    
  }
  
  # Show just the scatter plot
  else {
    
    return(list(dose_response_curve, NULL))
    
  }
}