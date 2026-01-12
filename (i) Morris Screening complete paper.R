library(deSolve) #version 1.40
library(ggplot2) #version 3.5.1
library(ggtext) #version 0.1.2
library(reshape2) #version 1.4.4


# parameters in the code are named differently than in paper, compare to Lever et al. 2016 cited in paper
# for updated parameter names, please refer to "new_parameter_names"


set.seed(250724)

#### Code part 1: Morris screening with evaluations ----
# the first part of the code runs the Morris screening for the defined options in the section "Morris Screening options"
# it will evaluate over the count defined and calculate confidence intervals

#### ODE System and parameters ----
ode_system <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    yYgen = 100
    yYdegr = 50000
    yPgen = 100
    yPdegr = 50000
    Pmax = 100
    
    dL <- - L * R * kon + koff * (C0 + C1 + C2)
    dR <- - L * R * kon + koff * (C0 + C1 + C2)
    dC0 <- kon * L * R - (koff + kp) * C0
    dC1 <- kp * C0 - (koff + phi * kp) * C1
    dC2 <- phi * kp * C1 - koff * C2
    dY <- yYgen * (Ymax - Y) - yYdegr * Y + lambda * C1 * (Ymax - Y)
    dP <- yPgen * (Pmax - P) - yPdegr * P + delta * Y * (Pmax - P) - my * C1 * P
    return(list(c(dL, dR, dC0, dC1, dC2, dY, dP)))
  })
}

parameter_bounds <- list(kon = c(0.0001, 100),   #copied koff Lever
                         koff = c(0.0001, 100),  #Lever
                         kp = c(0.00001, 1000),  #Lever
                         phi = c(0.00001, 1000), #Lever
                         Ymax = c(0.001, 1000),  #-3 to 3 to cover log-space
                         lambda = c(1, 1000000), #Lever
                         delta = c(1, 1000000),  #Lever
                         my = c(1, 1000000))     #Lever

#### Morris Screening options----
p <- 4  # number of levels in p-level grid
delta_morris <- p / (2 * (p - 1))  # step size delta
num_trajectories <- 80  # number of trajectories
Evaluations <- 10 # number of evaluations

#### generating trajectories ----
# Generate trajectory based on the specified sampling method
generate_trajectory <- function(parameter_bounds, p, delta_morris) {
  # Define the levels in the log space
  levels <- lapply(parameter_bounds, function(bounds) {
    log_bounds <- log(bounds)
    c(exp(log_bounds[1]), 
      exp(log_bounds[1] + (log_bounds[2] - log_bounds[1]) / 3), 
      exp(log_bounds[1] + 2 * (log_bounds[2] - log_bounds[1]) / 3), 
      exp(log_bounds[2]))
  })
  
  # Sample initial values
  sampled_values <- sapply(levels, function(level) sample(level, 1))
  
  # Initialize trajectory matrix
  trajectory <- matrix(NA, nrow = length(sampled_values) + 1, ncol = length(sampled_values))
  colnames(trajectory) <- names(parameter_bounds)
  trajectory[1, ] <- sampled_values
  
  # Adjust values based on delta_morris
  for (i in 1:length(sampled_values)) {
    new_point <- trajectory[i, ]
    if (new_point[i] == levels[[i]][1] || new_point[i] == levels[[i]][2]) {
      new_point[i] <- new_point[i] * exp(delta_morris * (log(parameter_bounds[[i]][2]) - log(parameter_bounds[[i]][1])))
    } else {
      new_point[i] <- new_point[i] / exp(delta_morris * (log(parameter_bounds[[i]][2]) - log(parameter_bounds[[i]][1])))
    }
    trajectory[i + 1, ] <- new_point
  }
  
  return(trajectory)
}

# Function to find the stable value of P
find_stable_P <- function(P_values, tolerance = 1e-4, consecutive_points = 5) {
  for (i in 1:(length(P_values) - consecutive_points + 1)) {
    if (max(abs(diff(P_values[i:(i + consecutive_points - 1)]))) < tolerance) {
      return(P_values[i + consecutive_points - 1])
    }
  }
  return(NA)  # Return NA if no stable value is found
}

# Function to calculate the elementary effects (EE)
calculate_EE <- function(trajectory) {
  EE <- list()
  
  for (i in 1:(nrow(trajectory) - 1)) {
    params_1 <- as.list(trajectory[i, ])
    params_2 <- as.list(trajectory[i + 1, ])
    
    out_1 <- ode(y = initial_state, times = times, func = ode_system, parms = params_1)
    out_2 <- ode(y = initial_state, times = times, func = ode_system, parms = params_2)
    
    P_1 <- out_1[, "P"]
    P_2 <- out_2[, "P"]
    
    stable_P_1 <- find_stable_P(P_1)
    stable_P_2 <- find_stable_P(P_2)
    
    if (!is.na(stable_P_1) && !is.na(stable_P_2)) {
      EE[[colnames(trajectory)[i]]] <- abs(stable_P_2 - stable_P_1)
    } else {
      EE[[colnames(trajectory)[i]]] <- NA  # If no stable value is found, assign NA
    }
  }
  
  return(EE)
}

# Custom function to calculate mean and standard deviation
calculate_stats <- function(ee_values) {
  # Remove NA values
  ee_values <- ee_values[!is.na(ee_values)]
  r <- length(ee_values)
  
  if (r == 0) {
    return(list(mean = NA, sd = NA))
  }
  
  # Calculate mean
  mu <- sum(abs(ee_values)) / r
  
  # Calculate standard deviation
  sigma <- sqrt(sum((abs(ee_values) - mu)^2) / (r - 1))
  
  return(list(mean = mu, sd = sigma))
}

# Define initial state and times for ODE solving
initial_state <- c(L = 50, R = 50, C0 = 0, C1 = 0, C2 = 0, Y = 0, P = 0)
times <- seq(0, 5000, by = 1) #simulate until 5000 so >99% of parameter sets run into steady state

# Repeat the entire process for the number of evaluations
all_EE_stats <- list()

for (eval in 1:Evaluations) {
  # Generate trajectories with different start values
  trajectories <- lapply(1:num_trajectories, function(x) generate_trajectory(parameter_bounds, p, delta_morris))
  
  # Calculate EE for each trajectory
  EE_results <- lapply(trajectories, calculate_EE)
  
  # Combine results for each parameter
  EE_combined <- list()
  for (param in names(parameter_bounds)) {
    EE_combined[[param]] <- sapply(EE_results, function(res) res[[param]])
  }
  
  
  # Calculate mean and standard deviation for each parameter
  EE_stats <- lapply(EE_combined, calculate_stats)
  
  # Store the results for this evaluation, including the raw EE values
  all_EE_stats[[eval]] <- list(EE_combined = EE_combined, EE_stats = EE_stats)
}

# Output the EE results with mean and standard deviation for each evaluation
all_EE_stats


#### calculate 95% confidence intervals ----
compute_ci <- function(values, alpha = 0.05) {
  vals <- values[!is.na(values)]
  n <- length(vals)
  m <- if (n > 0) mean(vals) else NA_real_ #calculate mean
  if (n < 2) {
    return(list(mean = m, ci_low = NA_real_, ci_high = NA_real_,
                halfwidth = NA_real_, rel_halfwidth = NA_real_, n = n))
  }
  s <- sd(vals) #standard deviation calc
  se <- s / sqrt(n) #standard error calc
  tcrit <- qt(1 - alpha/2, df = n - 1) # calculate critical t value
  hw <- tcrit * se # calculte halfwidth
  rel <- if (!is.na(m) && m != 0) hw / abs(m) else NA_real_ #calculate relative halfwidth
  list(mean = m, ci_low = m - hw, ci_high = m + hw,
       halfwidth = hw, rel_halfwidth = rel, n = n)
}

#create df with parameter + mean, sd and ci
make_plot_data_ci <- function(all_EE_stats, alpha = 0.05) {
  plot_data <- data.frame(Parameter = character(),
                          Mu = numeric(), Sigma = numeric(),
                          Mu_low = numeric(), Mu_high = numeric(),
                          Sigma_low = numeric(), Sigma_high = numeric(),
                          stringsAsFactors = FALSE)
  for (param in names(parameter_bounds)) {
    mu_values <- sapply(all_EE_stats, function(x) x$EE_stats[[param]]$mean)
    sigma_values <- sapply(all_EE_stats, function(x) x$EE_stats[[param]]$sd)
    ci_mu <- compute_ci(mu_values, alpha)
    ci_sigma <- compute_ci(sigma_values, alpha)
    plot_data <- rbind(plot_data, data.frame(
      Parameter = param,
      Mu = ci_mu$mean, Sigma = ci_sigma$mean,
      Mu_low = ci_mu$ci_low, Mu_high = ci_mu$ci_high,
      Sigma_low = ci_sigma$ci_low, Sigma_high = ci_sigma$ci_high
    ))
  }
  plot_data
}

ci_summary_table <- function(all_EE_stats, alpha = 0.05) {
  out <- data.frame(Parameter = character(),
                    Mu_mean = numeric(), Mu_CI_low = numeric(), Mu_CI_high = numeric(), Mu_rel_halfwidth = numeric(),
                    Sigma_mean = numeric(), Sigma_CI_low = numeric(), Sigma_CI_high = numeric(), Sigma_rel_halfwidth = numeric(),
                    n = integer(),
                    stringsAsFactors = FALSE)
  for (param in names(parameter_bounds)) {
    mu_values <- sapply(all_EE_stats, function(x) x$EE_stats[[param]]$mean)
    sigma_values <- sapply(all_EE_stats, function(x) x$EE_stats[[param]]$sd)
    ci_mu <- compute_ci(mu_values, alpha)
    ci_sigma <- compute_ci(sigma_values, alpha)
    out <- rbind(out, data.frame(
      Parameter = param,
      Mu_mean = ci_mu$mean, Mu_CI_low = ci_mu$ci_low, Mu_CI_high = ci_mu$ci_high, Mu_rel_halfwidth = ci_mu$rel_halfwidth,
      Sigma_mean = ci_sigma$mean, Sigma_CI_low = ci_sigma$ci_low, Sigma_CI_high = ci_sigma$ci_high, Sigma_rel_halfwidth = ci_sigma$rel_halfwidth,
      n = ci_mu$n
    ))
  }
  out
}




#### Plotting ----

#### create ci (confidence interval) data and plot ----
plot_data_ci <- make_plot_data_ci(all_EE_stats, alpha = 0.05)

ggplot(plot_data_ci, aes(x = Mu, y = Sigma, label = Parameter)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_text(vjust = -1) +
  geom_errorbar(aes(ymin = Sigma_low, ymax = Sigma_high), width = 0.2, color = "black") +
  geom_errorbarh(aes(xmin = Mu_low, xmax = Mu_high), height = 0.2, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    title = "Morris Screening with 95% confidence intervals",
    x = expression(mu*"*"),
    y = expression(sigma)
  ) +
  xlim(0, NA) +
  ylim(0, NA)



#plot for paper 
ggplot(plot_data_ci, aes(x = Mu, y = Sigma, label = "")) +
  geom_point(size = 3, color = "deepskyblue2") +
  geom_text(vjust = -1) +
  geom_errorbar(aes(ymin = Sigma_low, ymax = Sigma_high), width = 0.2, color = "black") +
  geom_errorbarh(aes(xmin = Mu_low, xmax = Mu_high), height = 0.2, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    x = expression(mu*"* (absolute mean of elementary effects)"),
    y = expression(sigma*" (standard deviation of elementary effects)")) +
  scale_y_continuous(
    limits = c(0, 40),
    breaks = seq(0, 40, by = 10),
    expand = expansion(mult = 0, add = 0)
  )+
  scale_x_continuous(
    limits = c(0, 40),
    breaks = seq(0, 40, by = 10),
    expand = expansion(mult = 0, add = 0)
  )+
  theme(
    axis.text.x = element_markdown(size = 16, margin = margin(t = 5)),
    axis.text.y = element_markdown(size = 16, margin = margin(r = 5)),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.margin = margin(10, 20, 10, 20)
  )



# ci summary
summary_ci <- ci_summary_table(all_EE_stats, alpha = 0.05)
print(summary_ci)


# calculate CIs based on cumulative evaluations -> one evaluation after another is added and hw is saved
ci_changes <- data.frame(Evaluations = integer(), Parameter = character(),
                         Mu_halfwidth = numeric(), Sigma_halfwidth = numeric())

for (num_evals in 2:Evaluations) {
  
  cumulative_stats <- all_EE_stats[1:num_evals]
  
  
  ci_summary <- ci_summary_table(cumulative_stats, alpha = 0.05)
  
  
  for (param in ci_summary$Parameter) {
    mu_hw <- ci_summary$Mu_rel_halfwidth[ci_summary$Parameter == param]
    sigma_hw <- ci_summary$Sigma_rel_halfwidth[ci_summary$Parameter == param]
    
    ci_changes <- rbind(ci_changes, data.frame(Evaluations = num_evals, Parameter = param,
                                               Mu_halfwidth = mu_hw, Sigma_halfwidth = sigma_hw))
  }
}

new_parameter_names <- c(  "kon" = "k<sub>on</sub>",
                           "koff" = "k<sub>off</sub>",
                           "kp" = "k<sub>p</sub>",
                           "phi" = "k<sub>lim</sub>",
                           "Ymax" = "M<sub>max</sub>",
                           "lambda" = "k<sub>amp</sub>",
                           "delta" = "k<sub>act</sub>",
                           "my" = "k<sub>inh</sub>")

ci_changes$Parameter <- factor(
  ci_changes$Parameter,
  levels = names(new_parameter_names),
  labels = unname(new_parameter_names)
)


# plot relative halfwidth over number of evaluations
ggplot(ci_changes, aes(x = Evaluations)) +
  geom_line(linewidth = 1, aes(y = Mu_halfwidth, color = Parameter, linetype = "Mu halfwidth")) +
  geom_line(linewidth = 1, aes(y = Sigma_halfwidth, color = Parameter, linetype = "Sigma halfwidth")) +
  theme_minimal() +
  labs(
    # title = "Change of confidence interval halfwidth over evaluations",
    x = "Number of evaluations",
    y = "Relative halfwidth of confidence intervals",
    linetype = "Type"
  ) +
  scale_linetype_manual(values = c("Mu halfwidth" = "solid", "Sigma halfwidth" = "dashed"),
                        breaks = c("Mu halfwidth", "Sigma halfwidth"),
                        labels = c("\u03BC *", "\u03C3"))+
  theme(
    axis.text.y = element_markdown(size = 18, margin = margin(r = 8)),     
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.title = element_text(size = 14), 
    legend.text = element_markdown(size = 14),   
    legend.key.size = unit(1.5, "cm")   
  )



# check for hw changes below tolerance
check_stable_ci_changes <- function(ci_changes, tolerance = 0.05, consecutive_points = 3) {
  stable_params <- list()
  
  for (param in unique(ci_changes$Parameter)) {
    param_data <- ci_changes[ci_changes$Parameter == param, ]
    mu_halfwidths <- param_data$Mu_halfwidth
    sigma_halfwidths <- param_data$Sigma_halfwidth
    
    mu_stable <- FALSE
    sigma_stable <- FALSE
    mu_stable_at <- NA
    sigma_stable_at <- NA
    
    for (i in 1:(length(mu_halfwidths) - consecutive_points + 1)) {
      if (max(abs(diff(mu_halfwidths[i:(i + consecutive_points - 1)]))) < tolerance) {
        mu_stable <- TRUE
        mu_stable_at <- param_data$Evaluations[i + consecutive_points - 1]
        break
      }
    }
    
    for (i in 1:(length(sigma_halfwidths) - consecutive_points + 1)) {
      if (max(abs(diff(sigma_halfwidths[i:(i + consecutive_points - 1)]))) < tolerance) {
        sigma_stable <- TRUE
        sigma_stable_at <- param_data$Evaluations[i + consecutive_points - 1]
        break
      }
    }
    
    if (mu_stable && sigma_stable) {
      stable_params[[param]] <- list(mu_stable_at = mu_stable_at, sigma_stable_at = sigma_stable_at)
    }
  }
  
  return(stable_params)
}

# Check: which parameter has a stable hw and when
stable_parameters <- check_stable_ci_changes(ci_changes, tolerance = 0.05, consecutive_points = 3)

if (length(stable_parameters) > 0) {
  for (param in names(stable_parameters)) {
    cat(sprintf("Parameter '%s' has a stable relative halfwidth at evaluation %d (Mu) und %d (Sigma).\n",
                param, stable_parameters[[param]]$mu_stable_at, stable_parameters[[param]]$sigma_stable_at))
  }
} else {
  print("No parameters with stable relative halfwidth")
}



#### Code part 2: generate evaluation for number of trajectories ----
# This part calculates the Morris indices for multiple trajectory counts 
# It keeps the old trajectories and adds 10 more each step to see how the sensitivity analysis stabilize 

set.seed(250724)
all_trajectories <- list()

# Initialize data frame to store results
results <- data.frame(Trajectories = integer(), Parameter = character(), Mu = numeric(), stringsAsFactors = FALSE)


# Loop to increase trajectories in steps of 10
for (num_trajectories in seq(10, 150, by = 10)) {
  # Generate new trajectories
  new_trajectories <- lapply(1:10, function(x) generate_trajectory(parameter_bounds, p, delta_morris))
  
  # Append new trajectories to the list of all trajectories
  all_trajectories <- c(all_trajectories, new_trajectories)
  
  # Calculate EE for each trajectory
  EE_results <- lapply(all_trajectories, calculate_EE)
  
  # Combine results for each parameter
  EE_combined <- list()
  for (param in names(parameter_bounds)) {
    EE_combined[[param]] <- sapply(EE_results, function(res) res[[param]])
  }
  
  # Calculate mean and standard deviation for each parameter
  EE_stats <- lapply(EE_combined, calculate_stats)
  
  # Store results in data frame
  for (param in names(EE_stats)) {
    results <- rbind(results, data.frame(Trajectories = length(all_trajectories), Parameter = param, Mu = EE_stats[[param]]$mean))
  }
}



results$Parameter <- factor(
  results$Parameter,
  levels = names(new_parameter_names),
  labels = unname(new_parameter_names)
)


# Plotting the results
ggplot(results, aes(x = Trajectories, y = Mu, color = Parameter, group = Parameter)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    #title = "Mu Values for Parameters with Increasing Number of Trajectories",
    x = "Number of Trajectories",
    y = expression(mu*"* (absolute mean)")
  ) +
  scale_color_discrete(name = "Parameter") +
  theme(
    axis.text.y = element_markdown(size = 18, margin = margin(r = 8)),     
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.title = element_text(size = 14), 
    legend.text = element_markdown(size = 14),   
    legend.key.size = unit(1.5, "cm")   
  )



# Calculate rankings
results$Rank <- ave(results$Mu, results$Trajectories, FUN = function(x) rank(-x, ties.method = "first"))

# Plotting the rankings
ggplot(results, aes(x = Trajectories, y = Rank, color = Parameter, group = Parameter)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Relative Ranking of Parameters with Increasing Number of Trajectories",
    x = "Number of Trajectories",
    y = "Rank"
  ) +
  scale_y_reverse(breaks = 1:length(parameter_bounds)) +
  scale_color_discrete(name = "Parameter") +
  theme(
    axis.text.y = element_markdown(size = 14, margin = margin(r = 8)),     
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.text = element_markdown(size = 14),   
    legend.key.size = unit(1, "cm")   
  )

# check for evalutions count with change in mu* below 20%
# returns the middle trajectory of the first stable window
check_stability_relative_extended_middle <- function(results, relative_threshold = 0.2, consecutive_points = 3, explain = TRUE) {
  unique_trajectories <- unique(results$Trajectories)
  
  for (i in seq_along(unique_trajectories)[-c(1:(consecutive_points - 1))]) {
    mu_values <- sapply(unique_trajectories[(i - consecutive_points + 1):i], function(traj) {
      results[results$Trajectories == traj, "Mu"]
    }, simplify = "array")
    
    # Relative changes for each parameter across the window
    relative_changes <- apply(mu_values, 1, function(mu) {
      abs(diff(mu)) / abs(mu[-length(mu)])
    })
    
    # All parameters must be below the threshold in all transitions of the window
    if (all(apply(relative_changes, 2, function(x) all(x < relative_threshold)))) {
      mid_idx <- i - floor((consecutive_points - 1) / 2)
      mid_traj <- unique_trajectories[mid_idx]
      if (explain) {
        message(sprintf("Returning the middle trajectory (%d) because it is the minimal n that confirms stability within the first stable window of 3 stable trajectories.", mid_traj))
      }
      return(mid_traj)
    }
  }
  return(NA)
}

stable_trajectory_number_relative_extended_middle <- check_stability_relative_extended_middle(results)
stable_trajectory_number_relative_extended_middle

