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


parameter_bounds <- list(kon = c(0.0001, 100),          #copied koff Lever
                         koff = c(0.0001, 100),         #Lever
                         kp = c(0.00001, 1000),         #Lever
                         phi = c(0.00001, 1000),        #Lever
                         lambda = c(1, 1000000),        #Lever
                         my = c(1, 1000000),            #Lever
                         deltaYmax = c(0.001, 1e9))     #lumped k_act * M_max (product of former bounds)

# fixed receptor amount; ligand axis 1e-4 to 1e4 covers L/R from 2e-6 to 200
R0_fixed <- 50

# original bounds of the lumped components k_act and M_max
# their product must equal parameter_bounds$deltaYmax
lumped_component_bounds <- list(k_act = c(1, 1000000),
                                M_max = c(0.001, 1000))

stopifnot(!anyNA(unlist(lumped_component_bounds)))
stopifnot(isTRUE(all.equal(
  lumped_component_bounds$k_act * lumped_component_bounds$M_max,
  parameter_bounds$deltaYmax
)))

# new parameter names (display labels, markdown)
new_parameter_names <- c("kon"       = "k<sub>on</sub>",
                         "koff"      = "k<sub>off</sub>",
                         "kp"        = "k<sub>p</sub>",
                         "phi"       = "k<sub>lim</sub>",
                         "lambda"    = "k<sub>amp</sub>",
                         "my"        = "k<sub>inh</sub>",
                         "deltaYmax" = "k<sub>act</sub>M<sub>max</sub>")

# shared colour palette for Sobol and Morris plots (identical block in both scripts)
# keyed by display label, so the same parameter always gets the same colour
param_colors <- c("kon"       = "#F8766D",
                  "koff"      = "#C49A00",
                  "kp"        = "#53B400",
                  "phi"       = "#00C094",
                  "lambda"    = "#00B6EB",
                  "my"        = "#A58AFF",
                  "deltaYmax" = "#FB61D7")
pal <- setNames(unname(param_colors[names(new_parameter_names)]),
                unname(new_parameter_names))
stopifnot(!anyNA(pal))


stopifnot(setequal(names(new_parameter_names), names(parameter_bounds)))

#### Morris Screening options----
p <- 64  # number of levels in p-level grid
delta_morris <- p / (2 * (p - 1))  # step size delta (in unit space, delta = 32/63 for p = 64)
num_trajectories <- 150  # number of trajectories
Evaluations <- 10 # number of evaluations

#### sampling transformation (identical to the Sobol script) ----
# inverse CDF of k_act * M_max with log-uniform k_act and M_max
# log10 of the product is the sum of two uniforms (trapezoidal distribution)
qlog_sum_uniform <- function(u, bounds1, bounds2) {
  lo <- log10(bounds1[1]) + log10(bounds2[1])
  w  <- sort(c(diff(log10(bounds1)), diff(log10(bounds2))))
  w1 <- w[1]
  w2 <- w[2]
  u1 <- w1 / (2 * w2)
  z  <- ifelse(u <= u1,
               sqrt(2 * w1 * w2 * u),
               ifelse(u <= 1 - u1,
                      u * w2 + w1 / 2,
                      w1 + w2 - sqrt(2 * w1 * w2 * pmax(1 - u, 0))))
  10^(lo + z)
}

# scale unit-cube values to log-uniform bounds
# deltaYmax is drawn as the product of log-uniform k_act and M_max
scale_to_bounds <- function(sample, bounds) {
  out <- mapply(function(x, range) 10^(log10(range[1]) + x * (log10(range[2]) - log10(range[1]))),
                as.data.frame(sample), bounds)
  out <- matrix(out, ncol = length(bounds))
  j <- which(names(bounds) == "deltaYmax")
  out[, j] <- qlog_sum_uniform(sample[, j],
                               lumped_component_bounds$k_act,
                               lumped_component_bounds$M_max)
  colnames(out) <- names(bounds)
  out
}

#### generating trajectories ----
# trajectories are built on the p-level grid {0, 1/(p-1), ..., 1} in the unit cube
# and then mapped to parameter values with scale_to_bounds()
# for log-uniform parameters this reproduces the equidistant log grid,
# for deltaYmax the levels are the corresponding quantiles of the product distribution
generate_trajectory <- function(parameter_bounds, p, delta_morris) {
  k <- length(parameter_bounds)
  step_levels <- round(delta_morris * (p - 1))   
  
  # sample initial levels (integer level indices avoid floating-point comparisons)
  level_idx <- matrix(NA_integer_, nrow = k + 1, ncol = k)
  level_idx[1, ] <- sample(0:(p - 1), k, replace = TRUE)
  
  # change one parameter per step: up if the step stays inside the grid, otherwise down
  for (i in 1:k) {
    new_point <- level_idx[i, ]
    if (new_point[i] + step_levels <= p - 1) {
      new_point[i] <- new_point[i] + step_levels
    } else {
      new_point[i] <- new_point[i] - step_levels
    }
    level_idx[i + 1, ] <- new_point
  }
  
  unit_trajectory <- level_idx / (p - 1)
  trajectory <- scale_to_bounds(unit_trajectory, parameter_bounds)
  colnames(trajectory) <- names(parameter_bounds)
  return(trajectory)
}

#### readout (identical to the Sobol script) ----
# analytical steady state of the Lever model for arbitrary L0 and R0
# total complex C = C0 + C1 + C2 solves C = K * (L0 - C) * (R0 - C); disc = discriminant, s = linear coefficient
# vectorized over L0
steady_state_P <- function(params, L0, R0 = R0_fixed) {
  with(as.list(params), {
    yYgen <- 100; yYdegr <- 50000; yPgen <- 100; yPdegr <- 50000; Pmax <- 100
    a  <- kon / (koff + kp)
    b  <- kp / (koff + phi * kp)
    c2 <- phi * kp / koff
    K  <- a * (1 + b + b * c2)
    s    <- K * (L0 + R0) + 1
    disc <- (K * (L0 - R0))^2 + 2 * K * (L0 + R0) + 1   # = s^2 - 4*K^2*L0*R0, always > 0
    C    <- 2 * K * L0 * R0 / (s + sqrt(disc))
    C1 <- b * C / (1 + b + b * c2)
    m  <- (yYgen + lambda * C1) / (yYgen + yYdegr + lambda * C1)   # m = M / M_max
    P  <- (yPgen + deltaYmax * m) * Pmax / (yPgen + yPdegr + deltaYmax * m + my * C1)
    unname(P)
  })
}

# maximum of the steady-state dose-response curve over L in [1e-4, 1e4] at R0 = R0_fixed
# coarse grid (as in the Shiny app) followed by local refinement with optimize()
logL_grid <- seq(-4, 4, by = 0.1)
max_P_over_L <- function(params, R0 = R0_fixed) {
  P_grid <- steady_state_P(params, L0 = 10^logL_grid, R0 = R0)
  j <- which.max(P_grid)
  if (j == 1 || j == length(logL_grid)) {
    return(c(P_max = P_grid[j], logL_opt = logL_grid[j], at_bound = 1))
  }
  opt <- optimize(function(x) steady_state_P(params, L0 = 10^x, R0 = R0),
                  interval = logL_grid[c(j - 1, j + 1)], maximum = TRUE, tol = 1e-8)
  c(P_max = opt$objective, logL_opt = opt$maximum, at_bound = 0)
}

# GSA readout: maximal steady-state response over the ligand axis
evaluate_readout <- function(params) {
  max_P_over_L(params)[["P_max"]]
}

# Function to calculate the elementary effects (EE)
# EE = |P(x + delta e_i) - P(x)| (absolute difference, not divided by delta)
calculate_EE <- function(trajectory) {
  # readout is evaluated once per trajectory point (k + 1 evaluations)
  P_traj <- apply(trajectory, 1, function(row) evaluate_readout(as.list(row)))
  stopifnot(all(is.finite(P_traj)))   # analytical solution, no non-finite readouts expected
  
  EE <- list()
  for (i in 1:(nrow(trajectory) - 1)) {
    EE[[colnames(trajectory)[i]]] <- abs(P_traj[i + 1] - P_traj[i])
  }
  return(EE)
}

# Custom function to calculate mean and standard deviation
calculate_stats <- function(ee_values) {
  stopifnot(length(ee_values) > 1, all(is.finite(ee_values)))
  r <- length(ee_values)
  
  # Calculate mean
  mu <- sum(abs(ee_values)) / r
  
  # Calculate standard deviation
  sigma <- sqrt(sum((abs(ee_values) - mu)^2) / (r - 1))
  
  return(list(mean = mu, sd = sigma))
}

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

# sanity check: all elementary effects finite (analytical solution)
n_EE_total <- sum(sapply(all_EE_stats, function(x) length(unlist(x$EE_combined))))
cat(sprintf("Valid elementary effects: %d of %d (100%% expected)\n",
            sum(sapply(all_EE_stats, function(x) sum(is.finite(unlist(x$EE_combined))))),
            n_EE_total))

# check: grid levels per parameter (log-uniform parameters = equidistant log grid)
grid_levels <- scale_to_bounds(matrix(rep((0:(p - 1)) / (p - 1), length(parameter_bounds)),
                                      ncol = length(parameter_bounds)),
                               parameter_bounds)
print(signif(grid_levels, 3))

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
plot_data_ci$Label <- unname(new_parameter_names[plot_data_ci$Parameter])

# common axis range for mu* and sigma (same scale on both axes so that the 1:1 line is meaningful)
axis_max    <- max(c(plot_data_ci$Mu_high, plot_data_ci$Sigma_high), na.rm = TRUE) * 1.05
axis_breaks <- pretty(c(0, axis_max))
axis_max    <- max(axis_breaks)

ggplot(plot_data_ci, aes(x = Mu, y = Sigma)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_richtext(aes(label = Label), vjust = -0.6, fill = NA, label.color = NA) +
  geom_errorbar(aes(ymin = Sigma_low, ymax = Sigma_high), width = 0.2, color = "black") +
  geom_errorbarh(aes(xmin = Mu_low, xmax = Mu_high), height = 0.2, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    title = "Morris Screening with 95% confidence intervals",
    subtitle = "dashed line: \u03C3/\u03BC* = 1",
    x = expression(mu*"*"),
    y = expression(sigma)
  ) +
  coord_cartesian(xlim = c(0, axis_max), ylim = c(0, axis_max))

#plot for paper
ggplot(plot_data_ci, aes(x = Mu, y = Sigma)) +
  geom_point(size = 3, color = "deepskyblue2") +
  geom_errorbar(aes(ymin = Sigma_low, ymax = Sigma_high), width = 0.2, color = "black") +
  geom_errorbarh(aes(xmin = Mu_low, xmax = Mu_high), height = 0.2, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    x = expression(mu*"* (absolute mean of elementary effects)"),
    y = expression(sigma*" (standard deviation of elementary effects)")) +
  scale_y_continuous(
    breaks = axis_breaks,
    expand = expansion(mult = 0, add = 0)
  )+
  scale_x_continuous(
    breaks = axis_breaks,
    expand = expansion(mult = 0, add = 0)
  )+
  coord_cartesian(xlim = c(0, axis_max), ylim = c(0, axis_max)) +
  theme(
    axis.text.x = element_markdown(size = 16, margin = margin(t = 5)),
    axis.text.y = element_markdown(size = 16, margin = margin(r = 5)),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.margin = margin(10, 20, 10, 20)
  )

# ci summary (printed with 3 significant digits)
summary_ci <- ci_summary_table(all_EE_stats, alpha = 0.05)
summary_ci$Mu_halfwidth    <- summary_ci$Mu_CI_high - summary_ci$Mu_mean
summary_ci$Sigma_halfwidth <- summary_ci$Sigma_CI_high - summary_ci$Sigma_mean
summary_ci_print <- summary_ci
num_cols <- vapply(summary_ci_print, is.numeric, logical(1)) & names(summary_ci_print) != "n"
summary_ci_print[num_cols] <- lapply(summary_ci_print[num_cols], signif, digits = 3)
print(summary_ci_print, row.names = FALSE)

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
  scale_color_manual(values = pal, name = "Parameter") +
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

# Returns the middle evaluation point (and optional value) of the first stable window
check_stable_ci_changes_middle <- function(ci_changes, tolerance = 0.05, consecutive_points = 3) {
  stable_params <- list()
  
  for (param in unique(ci_changes$Parameter)) {
    param_data <- ci_changes[ci_changes$Parameter == param, ]
    mu_halfwidths <- param_data$Mu_halfwidth
    sigma_halfwidths <- param_data$Sigma_halfwidth
    
    mu_stable_at <- NA
    sigma_stable_at <- NA
    mu_mid_value <- NA
    sigma_mid_value <- NA
    
    n <- nrow(param_data)
    win_len <- consecutive_points
    if (n >= win_len) {
      for (i in 1:(n - win_len + 1)) {
        mu_window <- mu_halfwidths[i:(i + win_len - 1)]
        sigma_window <- sigma_halfwidths[i:(i + win_len - 1)]
        
        mu_ok <- max(abs(diff(mu_window))) < tolerance
        sigma_ok <- max(abs(diff(sigma_window))) < tolerance
        
        # Middle index within the window
        mid_offset <- floor((win_len - 1) / 2)
        mid_idx <- i + mid_offset
        
        if (mu_ok && is.na(mu_stable_at)) {
          mu_stable_at <- param_data$Evaluations[mid_idx]
          mu_mid_value <- mu_halfwidths[mid_idx]
        }
        if (sigma_ok && is.na(sigma_stable_at)) {
          sigma_stable_at <- param_data$Evaluations[mid_idx]
          sigma_mid_value <- sigma_halfwidths[mid_idx]
        }
        
        if (!is.na(mu_stable_at) && !is.na(sigma_stable_at)) break
      }
    }
    
    if (!is.na(mu_stable_at) && !is.na(sigma_stable_at)) {
      stable_params[[param]] <- list(
        mu_stable_at = mu_stable_at,
        sigma_stable_at = sigma_stable_at,
        mu_mid_value = mu_mid_value,
        sigma_mid_value = sigma_mid_value
      )
    }
  }
  
  return(stable_params)
}

# returns the middle evaluation count of the first stable window
stable_parameters_middle <- check_stable_ci_changes_middle(ci_changes, tolerance = 0.05, consecutive_points = 3)

if (length(stable_parameters_middle) > 0) {
  for (param in names(stable_parameters_middle)) {
    cat(sprintf(
      "Parameter '%s' has stable relative halfwidths at the middle point: evaluation %d (Mu, value=%.6f) and %d (Sigma, value=%.6f).\n",
      param,
      stable_parameters_middle[[param]]$mu_stable_at, stable_parameters_middle[[param]]$mu_mid_value,
      stable_parameters_middle[[param]]$sigma_stable_at, stable_parameters_middle[[param]]$sigma_mid_value
    ))
  }
} else {
  print("No parameters with stable relative halfwidth")
}



#### Code part 2: Morris sample size evaluation ----
# mu* is calculated for increasing numbers of trajectories (old trajectories are kept, new ones are added)
# convergence: relative change of the mu* below 20 % over the
# final window of 3 trajectory counts, checked with the same function as in the Sobol script
# (check_final_stability, rel_tolerance = 0.2, abs_tolerance = 0, window_size = 3)
# NOTE: requires generate_trajectory(), calculate_EE() and new_parameter_names from code part 1

set.seed(250724)


trajectory_steps <- seq(10, 170, by = 10)

EE_matrix <- matrix(NA_real_, nrow = 0, ncol = length(parameter_bounds),
                    dimnames = list(NULL, names(parameter_bounds)))
results <- data.frame()

for (n_traj in trajectory_steps) {
  # add only the new trajectories; EE of earlier trajectories are kept
  n_new  <- n_traj - nrow(EE_matrix)
  new_EE <- t(vapply(seq_len(n_new), function(x) {
    unlist(calculate_EE(generate_trajectory(parameter_bounds, p, delta_morris)))[names(parameter_bounds)]
  }, numeric(length(parameter_bounds))))
  EE_matrix <- rbind(EE_matrix, new_EE)
  
  mu_star <- colMeans(abs(EE_matrix))
  sigma   <- apply(EE_matrix, 2, sd)
  results <- rbind(results, data.frame(Trajectories = n_traj,
                                       Parameter    = names(parameter_bounds),
                                       Mu           = mu_star,
                                       Sigma        = sigma,
                                       Mu_norm      = mu_star / sum(mu_star)))
}

# convergence check (same function as in the Sobol script, 20 % relative tolerance) 
rel_tolerance <- 0.2
abs_tolerance <- 0
window_size   <- 3

# Check whether the final window of sample sizes is stable
check_final_stability <- function(
    N,
    values,
    rel_tolerance = 0.01,
    abs_tolerance = 0.005,
    window_size = 5
) {
  
  stopifnot(length(N) == length(values))
  stopifnot(length(values) >= window_size)
  
  ord <- order(N)
  N <- N[ord]
  values <- values[ord]
  
  idx <- (length(values) - window_size + 1):length(values)
  reference <- values[idx[1]]
  
  allowed_deviation <- max(
    abs_tolerance,
    rel_tolerance * abs(reference)
  )
  
  absolute_deviation <- abs(values[idx] - reference)
  converged <- all(absolute_deviation <= allowed_deviation)
  
  data.frame(
    Converged = converged,
    Final_window_start_N = if (converged) N[idx[1]] else NA_integer_,
    Window_end_N = N[idx[length(idx)]],
    Reference = reference,
    Allowed_abs_deviation = allowed_deviation,
    Max_abs_deviation = max(absolute_deviation)
  )
}

stability_mu <- do.call(
  rbind,
  lapply(names(parameter_bounds), function(parameter) {
    d <- results[results$Parameter == parameter, ]
    cbind(
      Index = "mu*",
      Parameter = parameter,
      check_final_stability(
        N = d$Trajectories,
        values = d$Mu,
        rel_tolerance = rel_tolerance,
        abs_tolerance = abs_tolerance,
        window_size = window_size
      )
    )
  })
)

cat("\nAll mu* converged (relative change < 20 %):", all(stability_mu$Converged), "\n\n")
print(stability_mu[, c("Index", "Parameter", "Converged", "Final_window_start_N", "Window_end_N",
                       "Max_abs_deviation", "Allowed_abs_deviation")], row.names = FALSE)

# ranking at the largest number of trajectories (compare with the STi ranking of the Sobol script)
final_mu <- results[results$Trajectories == max(results$Trajectories), ]
print(final_mu[order(final_mu$Mu, decreasing = TRUE), c("Parameter", "Mu", "Mu_norm", "Sigma")], row.names = FALSE)

# legend order = descending mu* at the largest number of trajectories
legend_order <- as.character(final_mu$Parameter[order(final_mu$Mu, decreasing = TRUE)])
results$Parameter <- factor(
  results$Parameter,
  levels = legend_order,
  labels = unname(new_parameter_names[legend_order])
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
  scale_color_manual(values = pal, name = "Parameter") +
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
