
view <- function(data, x_var = 'time', y_var, visualization = 'single', start, end) {
  
  # #######################################################################################
  # Description: show the experimental data, color-coded according to the condition
  #
  # Inputs: 
  #         data:           data-frame obtained from 'cleaning_data.R' <type:tibble>
  #         parameter:      model parameter, which can be
  #                         'od':               optical density
  #                         'flu':              fluorescence value
  #                         'growth_rate':      growth rate
  #                         'phi':              heterologous fraction
  #                         'production_rate':  production rate
  #         visualization:  type of plot generated by the function
  #                         'single':           legend-based plot (default value)
  #                         'multiple':         facet_wrap-based plot
  #
  # Output: 
  #         compare_plot: comparison between Condition concentrations <type:ggplot2>
  #
  # #######################################################################################
  
  # Transform conditions to factors
  data <- data %>%
    mutate(Condition = as.factor(Condition))
  
  # Define units for each variable
  ylab_expression <- switch(
    y_var,
    od = expression(OD ~ (A)),
    flu = expression(FU ~ (A)),
    growth_rate = expression(µ ~ (h^-1)),
    phi = expression(Φ[H] ~ (FU) ~ (OD^-1)),
    production_rate = expression(ρ[H] ~ (h^-1) ~ (OD^-1)),
    default = expression("Default Label")  # Set a default label if no matches are found
  )
  
  # 'single' view: legend-based plot 
  if (visualization == 'single') {
    ggplot(data = data, aes_string(x = x_var, y = y_var, color = 'Condition')) +
      geom_point() +
      geom_smooth(se = FALSE, span = 0.2) +
      geom_vline(xintercept = c(start, end), linetype = "dashed", color = "black",  linewidth = 1) +
      xlab('time (h)') +
      ylab(ylab_expression) +
      ggtitle(unique(data$strain)) +
      theme(plot.title = element_text(size = 20, face = "bold"))
  }
  
  # 'multiple' view: facet-wrap plot
  else if (visualization == 'multiple') {
    ggplot(data = data, aes_string(x = x_var, y = y_var, color = 'Condition')) +
      geom_point() +
      geom_smooth(se = FALSE, span = 0.2) +
      geom_vline(xintercept = c(start, end), linetype = "dashed", color = "black", linewidth = 1) +
      theme(legend.position = "none") +
      facet_wrap(~Condition) +
      xlab('time (h)') +
      ylab(ylab_expression) +
      ggtitle(unique(data$strain)) +
      theme(plot.title = element_text(size = 20, face = "bold"))
  }
}

accumulation <- function(name, inductor, time_interval, method) {
  
  # ##################################################################################
  # GFP_accumulation function:
  # 
  # Description: approximation (area under curve) of how GFP molecules does
  #              a strain have accumulated
  #
  # Inputs: 
  #         name:     data-frame obtained from 'cleaning_data.R' <type:tibble>
  #         Condition:     inductor concentration for analysis <unit:µM>
  #         method:   classic -> phi value at maximum growth rate as upper bound
  #                   maximum -> maximum phi value as upper bound
  #
  # Output: 
  #         p1:       growth rate plot <type:ggplot2>
  #         p2:       phi plot with area under curve highlighted to the time point
  #                   where growth rate is max <type:ggplot2>
  #
  # ##################################################################################
  
  strain <- name %>% 
    filter(Condition == inductor & time >= time_interval[[1]] & time <= time_interval[[2]])
  
  p1 <- ggplot(data = strain, aes(x = time, y = growth_rate)) +
    geom_point(color = "grey") +
    geom_smooth(se = FALSE, span = 0.2, color = "grey") +
    theme(legend.position = "none") +
    xlab('time (h)') +
    ylab(expression(µ~(h^-1)))+
    ggtitle(paste0("Growth rate (", inductor, ")"))+
    theme(plot.title = element_text(size = 15, face = "bold"))
  
  if (method == "bam") {
    
    p2 <- ggplot(data = strain, aes(x = time, y = phi)) +
      geom_point(color = "grey") +
      geom_smooth(se = FALSE, span = 0.2, color = "grey") +
      theme(legend.position = "none") +
      geom_ribbon(data = strain[1:match(max(strain$phi), strain$phi), ],
                  aes(x = time, ymax = phi), ymin = 0, fill = "green", alpha = 0.5) +
      xlab('time (h)') +
      ylab(expression(Φ[H] ~ (FU) ~ (OD^-1)))+
      ggtitle(paste0("Protein accumulation (", inductor, ")"))+
      theme(plot.title = element_text(size = 15, face = "bold"))
    
  }
  
  else if (method == "hwa") {
    
    # Find the index corresponding to the maximum growth rate
    max_index <- which.max(strain$growth_rate)
    
    # Calculate the maximum phi value up to the point of maximum growth rate
    max_phi <- max(strain$phi[1:max_index], na.rm = TRUE)
    
    # Plot strain$time vs strain$phi with geom_ribbon
    p2 <- ggplot(data = strain, aes(x = time, y = phi)) +
        geom_point(color = "grey") +
        geom_smooth(se = FALSE, span = 0.2, color = "grey") +
        theme(legend.position = "none") +
        geom_ribbon(data = subset(strain, time <= time[max_index]),
                    aes(x = time, ymax = phi, ymin = 0), fill = "green", alpha = 0.5) +
        xlab('time (h)') +
        ylab(expression(Φ[H] ~ (FU) ~ (OD^-1))) +
        ggtitle(paste0("Protein accumulation")) +
        theme(plot.title = element_text(size = 15, face = "bold"))
    
  }
  
  return(list(p1, p2))
}

hill_prediction <- function(data, 
                            model_fit,
                            time_interval,
                            growth_interval
                            ){
  
  # Define x-axis
  growth <- seq(growth_interval[[1]], growth_interval[[2]], by=0.1)
  
  # Extract concentrations
  aux <- data %>% 
    filter(Condition > min(Condition))
  
  #concentrations <- round(aux$Condition)
  concentrations <- aux$Condition
  
  # Phi predicted from hill function 
  predicted_fi <- model_fit %>%
    filter(Condition %in% concentrations) %>%
    select(Condition, fi, error)
  
  # Creating prediction lines with the hill parameters
  GP <- expand.grid(fi = predicted_fi$fi, gr = growth) %>%
    mutate(pr = fi*gr)
  
  GP_lines <- predicted_fi %>% left_join(GP, by = "fi") %>%
    mutate(lower = pr - qnorm(0.975)*error*gr,
           upper = pr + qnorm(0.975)*error*gr) %>% select(-error)
  
  # Comparing data
  #data$Condition <- round(data$Condition)
  
  compare_plot <- data %>%
    filter(Condition > min(Condition) & time >= time_interval[[1]] & time <= time_interval[[2]]) %>%
    mutate(Condition = factor(Condition)) %>%
    ggplot() +
    geom_point(aes(growth_rate, production_rate/1000), color = "black", fill = "white", shape = 21, size = 2.5, stroke = 1) + # Increase size of geom_point
    geom_smooth(aes(growth_rate, production_rate/1000), method = "lm", se = TRUE, fill = "grey80", color = "black", linetype = "dashed") + # Black dashed geom_smooth
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_line(data = GP_lines %>% mutate(Condition = factor(Condition)), aes(gr, pr/1000), color = "black", size = 1) + # Black lines
    facet_wrap(~Condition) +
    xlab(expression(Growth~rate~(h^-1))) +
    ylab(expression(Production~rate~(10^3*FU~OD[600]^-1*h^-1))) +
    guides(color = FALSE)  # Remove the legend for color

  return(compare_plot)
  
}