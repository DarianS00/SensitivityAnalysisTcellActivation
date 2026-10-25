library(deSolve) #version 1.40
library(doParallel) #version 1.0.17
library(foreach) #version 1.5.2
library(ggplot2) #version 3.5.1
library(ggtext) #version 0.1.2
library(randtoolbox) #version 2.0.5
library(svglite) #version 2.2.2



# parameters in the code are named differently than in paper, compare to Lever et al. 2016 cited in paper
# for updated parameter names, please refer to "new_parameter_names"


set.seed(250724)


#### Code part 1: Sobol screening for first-order, total-order and total-interaction indices ----
# the first part of the code runs the Sobol screening for the defined sample size (N) 

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
    # Y = M / M_max (normalized mediator); M_max cancels, only deltaYmax = k_act * M_max remains
    dY <- yYgen * (1 - Y) - yYdegr * Y + lambda * C1 * (1 - Y)
    dP <- yPgen * (Pmax - P) - yPdegr * P + deltaYmax * Y * (Pmax - P) - my * C1 * P
    return(list(c(dL, dR, dC0, dC1, dC2, dY, dP)))
  })
}

parameter_bounds <- list(kon = c(0.0001, 100),          #copied koff Lever
                         koff = c(0.0001, 100),         #Lever
                         kp = c(0.00001, 1000),         #Lever
                         phi = c(0.00001, 1000),        #Lever
                         lambda = c(1, 1000000),        #Lever
                         my = c(1, 1000000),            #Lever
                         deltaYmax = c(0.001, 1e9))     #lumped k_act * M_max (product of former bounds)



# number of samples (N) and parameters (k)
N <- 2^16  
k <- length(parameter_bounds) 

# generate sobol sequence for matrix A and B (2 * k)
sobol_samples <- sobol(N, 2 * k)

# first k columns for A, second k columns for B
A <- sobol_samples[, 1:k]
B <- sobol_samples[, (k + 1):(2 * k)]

# original bounds of the lumped components k_act and M_max
# fill in the values; their product must equal parameter_bounds$deltaYmax
lumped_component_bounds <- list(k_act = c(1, 1000000),
                                M_max = c(0.001, 1000))

stopifnot(!anyNA(unlist(lumped_component_bounds)))
stopifnot(isTRUE(all.equal(
  lumped_component_bounds$k_act * lumped_component_bounds$M_max,
  parameter_bounds$deltaYmax
)))

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
                      w1 + w2 - sqrt(2 * w1 * w2 * (1 - u))))
  10^(lo + z)
}

# scale Sobol samples to log-uniform bounds
# deltaYmax is drawn as the product of log-uniform k_act and M_max
scale_to_bounds <- function(sample, bounds) {
  out <- mapply(function(x, range) 10^(log10(range[1]) + x * (log10(range[2]) - log10(range[1]))),
                as.data.frame(sample), bounds)
  j <- which(names(bounds) == "deltaYmax")
  out[, j] <- qlog_sum_uniform(sample[, j],
                               lumped_component_bounds$k_act,
                               lumped_component_bounds$M_max)
  out
}

A <- as.data.frame(scale_to_bounds(A, parameter_bounds))
B <- as.data.frame(scale_to_bounds(B, parameter_bounds))

colnames(A) <- colnames(B) <- names(parameter_bounds)  # Namen setzen


# generate lsit of ABi matrices (AB_list contains k matrices, where one column from A has been replaced with the corresp. column fromn B)
AB_list <- lapply(1:k, function(i) {
  AB_i <- A  # copy A
  AB_i[, i] <- B[, i]  # replace column i with values from B
  return(AB_i)
})

# analytical steady state of the Lever model for arbitrary L0 and R0
# total complex C = C0 + C1 + C2 solves C = K * (L0 - C) * (R0 - C); disc = discriminant, s = linear coefficient 
# vectorized over L0; constants must match ode_system
steady_state_P <- function(params, L0, R0 = 50) {
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

# maximum of the steady-state dose-response curve over L in [1e-4, 1e4], R0 = 50
# coarse grid (as in the Shiny app) followed by local refinement with optimize()
logL_grid <- seq(-4, 4, by = 0.1)
max_P_over_L <- function(params, R0 = 50) {
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

time1 <- Sys.time()


# number of cores for parallel calculation
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

worker_exports <- c(
  "evaluate_readout",
  "max_P_over_L",
  "steady_state_P",
  "logL_grid"
)

# parallel calculation for matrix A
f_A <- foreach(i = 1:N, .combine = c, .export = worker_exports) %dopar% {
  evaluate_readout(A[i, ])
}

# parallel calculation for matrix B
f_B <- foreach(i = 1:N, .combine = c, .export = worker_exports) %dopar% {
  evaluate_readout(B[i, ])
}

# parallel calculation for ABi matrices
f_AB <- foreach(i = 1:k, .combine = 'c', .multicombine = TRUE,  .export = worker_exports) %dopar% {
  list(apply(AB_list[[i]], 1, evaluate_readout))
}

stopCluster(cl)

time2 <- Sys.time()
time2-time1

# sanity check: all readouts finite (analytical solution, no convergence failures expected)
n_nonfinite <- sum(!is.finite(c(f_A, f_B, unlist(f_AB))))
cat(sprintf("Non-finite readouts: %d of %d\n", n_nonfinite, N * (k + 2)))


# calculate model variance 
f_all <- c(f_A, f_B)
V_Y <- mean(f_all^2, na.rm = TRUE) - mean(f_all, na.rm = TRUE)^2


#### calculate sensitivity indices ----

#first order index
Vi <- sapply(1:k, function(i) {
  mean(f_B * (f_AB[[i]] - f_A), na.rm = TRUE) 
})


Si <- sapply(1:k, function(i) {
  Vi[[i]] / V_Y
})


parameter_names <- names(parameter_bounds)

# total order index
STi <- sapply(1:k, function(i) {
  (mean((f_A - f_AB[[i]])^2, na.rm = TRUE) / 2) / V_Y
})


# total interaction index D_ij = STi + STj - ST{i,j}
# ST{i,j} = total effect of the group {i,j}, Jansen-type estimator from AB_i and AB_j
# D_ij = sum of all Sobol terms containing both i and j (upper bound for the pure S_ij)
STij <- matrix(NA_real_, nrow = k, ncol = k)
Dij  <- matrix(NA_real_, nrow = k, ncol = k)
for (i in 1:(k - 1)) {
  for (j in (i + 1):k) {
    STij[i, j] <- (mean((f_AB[[i]] - f_AB[[j]])^2) / 2) / V_Y
    Dij[i, j]  <- STi[i] + STi[j] - STij[i, j]
    STij[j, i] <- STij[i, j]
    Dij[j, i]  <- Dij[i, j]
  }
}
rownames(Dij) <- colnames(Dij) <- parameter_names




# new parameter names 
new_parameter_names <- c("kon"       = "k<sub>on</sub>",
                         "koff"      = "k<sub>off</sub>",
                         "kp"        = "k<sub>p</sub>",
                         "phi"       = "k<sub>lim</sub>",
                         "lambda"    = "k<sub>amp</sub>",
                         "deltaYmax" = "k<sub>act</sub>M<sub>max</sub>",
                         "my"        = "k<sub>inh</sub>")

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

# bootstrap confidence intervals for Si and STi 
# resample the N rows (same rows in A, B and all AB_i) with replacement
# and recompute V_Y, Si and STi from the existing model outputs
n_boot <- 1000
set.seed(250724)

Si_boot  <- matrix(NA_real_, nrow = n_boot, ncol = k)
STi_boot <- matrix(NA_real_, nrow = n_boot, ncol = k)
pairs    <- t(combn(k, 2))
Dij_boot <- matrix(NA_real_, nrow = n_boot, ncol = nrow(pairs))

for (b in 1:n_boot) {
  idx     <- sample(N, N, replace = TRUE)
  fA_b    <- f_A[idx]
  fB_b    <- f_B[idx]
  fAB_b   <- lapply(f_AB, function(x) x[idx])
  f_all_b <- c(fA_b, fB_b)
  V_Y_b   <- mean(f_all_b^2, na.rm = TRUE) - mean(f_all_b, na.rm = TRUE)^2
  
  Si_boot[b, ]  <- sapply(1:k, function(i) {
    mean(fB_b * (fAB_b[[i]] - fA_b), na.rm = TRUE) / V_Y_b
  })
  STi_boot[b, ] <- sapply(1:k, function(i) {
    (mean((fA_b - fAB_b[[i]])^2, na.rm = TRUE) / 2) / V_Y_b
  })
  Dij_boot[b, ] <- apply(pairs, 1, function(pr) {
    STij_b <- (mean((fAB_b[[pr[1]]] - fAB_b[[pr[2]]])^2) / 2) / V_Y_b
    STi_boot[b, pr[1]] + STi_boot[b, pr[2]] - STij_b
  })
}



# 95% percentile intervals (column 1 = lower, column 2 = upper)
Si_ci  <- t(apply(Si_boot,  2, quantile, probs = c(0.025, 0.975)))
STi_ci <- t(apply(STi_boot, 2, quantile, probs = c(0.025, 0.975)))

# table for the manuscript: estimate [lower, upper], sorted by STi
fmt_ci <- function(m, lo, hi) {
  sprintf(
    "%s [%s, %s]",
    formatC(m, digits = 3, format = "g"),
    formatC(lo, digits = 3, format = "g"),
    formatC(hi, digits = 3, format = "g")
  )
}

sobol_table <- data.frame(
  Parameter = parameter_names,
  Si  = fmt_ci(Si,  Si_ci[, 1],  Si_ci[, 2]),
  STi = fmt_ci(STi, STi_ci[, 1], STi_ci[, 2]),
  STi_value = STi
)
sobol_table <- sobol_table[order(sobol_table$STi_value, decreasing = TRUE), ]
sobol_table$STi_value <- NULL
print(sobol_table, row.names = FALSE)

Dij_ci <- t(apply(Dij_boot, 2, quantile, probs = c(0.025, 0.975)))
Dij_lo <- Dij_hi <- matrix(NA_real_, nrow = k, ncol = k, dimnames = dimnames(Dij))
Dij_lo[pairs] <- Dij_ci[, 1]; Dij_lo[pairs[, 2:1]] <- Dij_ci[, 1]
Dij_hi[pairs] <- Dij_ci[, 2]; Dij_hi[pairs[, 2:1]] <- Dij_ci[, 2]
Dij_est <- Dij[pairs]
dij_table <- data.frame(
  Pair = paste(parameter_names[pairs[, 1]], parameter_names[pairs[, 2]], sep = " x "),
  Dij  = fmt_ci(Dij_est, Dij_ci[, 1], Dij_ci[, 2]),
  Dij_value = Dij_est
)
dij_table <- dij_table[order(dij_table$Dij_value, decreasing = TRUE), ]
dij_table$Dij_value <- NULL
print(dij_table, row.names = FALSE)



# create df for sensitivity indices with new parameter names
result <- data.frame(Parameter = new_parameter_names[names(parameter_bounds)], Si = Si, STi = STi,
                     Si_low = Si_ci[, 1], Si_high = Si_ci[, 2],
                     STi_low = STi_ci[, 1], STi_high = STi_ci[, 2])



# labels for point estimates and 95% bootstrap confidence intervals
fmt3 <- function(x) formatC(x, digits = 3, format = "f")

result$Si_value_label  <- fmt3(result$Si)
result$Si_ci_label     <- paste0(
  "[", fmt3(result$Si_low), ", ", fmt3(result$Si_high), "]"
)

result$STi_value_label <- fmt3(result$STi)
result$STi_ci_label    <- paste0(
  "[", fmt3(result$STi_low), ", ", fmt3(result$STi_high), "]"
)

# position labels slightly to the right of the upper CI limit
result$Si_label_pos  <- result$Si_high + 0.025
result$STi_label_pos <- result$STi_high + 0.025

# sort df for Si
result_sorted <- result[order(result$Si, decreasing = TRUE), ]

# Si barplot
ggplot(result_sorted, aes(x = reorder(Parameter, Si), y = Si)) +
  geom_bar(stat = "identity", fill = "deepskyblue", color = "black") +
  geom_errorbar(aes(ymin = Si_low, ymax = Si_high), width = 0.3) +
  geom_text(aes(y = Si_label_pos, label = Si_value_label), hjust = 0, vjust = -0.15, size = 4.5) +
  geom_text(aes(y = Si_label_pos, label = Si_ci_label), hjust = 0, vjust = 1.35, size = 3.5)+
  labs( x = "", y = "S<sub>i</sub>") +
  theme_light() +
  coord_flip( ylim = c(0, 1.20),expand = FALSE, clip = "off")+
  theme(
    axis.text.y = element_markdown(size = 18, margin = margin(r = 8)),     
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_markdown(size = 18),
    plot.margin = margin(10, 140, 10, 10)
    ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = 0, add = 0)  
  )

#sort df for STi
result_sorted_STi <- result[order(result$STi, decreasing = TRUE), ]

# STi barplot
ggplot(result_sorted_STi, aes(x = reorder(Parameter, STi), y = STi)) +
  geom_bar(stat = "identity", fill = "deepskyblue4", color = "black") +
  geom_errorbar(aes(ymin = STi_low, ymax = STi_high), width = 0.3) +
  geom_text(aes(y = STi_label_pos, label = STi_value_label), hjust = 0, vjust = -0.15, size = 4.5) +
  geom_text(aes(y = STi_label_pos, label = STi_ci_label), hjust = 0, vjust = 1.35, size = 3.5)+
  labs( x = "", y = "S<sub>Ti</sub>") +
  theme_light() +
  coord_flip( ylim = c(0, 1.20),expand = FALSE, clip = "off")+
  theme(
    axis.text.y = element_markdown(size = 18, margin = margin(r = 8)),     
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_markdown(size = 18),
    plot.margin = margin(10, 140, 10, 10)
    )+
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = 0, add = 0) 
  )


# heatmap of total interaction indices D_ij
# axes sorted by STi (descending), lower triangle only, diagonal = not applicable (grey)
# Dij_df is rebuilt from the Dij matrix here, so the block can be re-run safely
n_par <- length(parameter_names)
ord   <- parameter_names[order(STi, decreasing = TRUE)]

Dij_df <- as.data.frame(as.table(Dij), stringsAsFactors = FALSE)
colnames(Dij_df) <- c("Var1", "Var2", "value")
Dij_df$lo <- as.vector(Dij_lo)
Dij_df$hi <- as.vector(Dij_hi)

Dij_df$Var1 <- factor(new_parameter_names[Dij_df$Var1], levels = new_parameter_names[ord])
Dij_df$Var2 <- factor(new_parameter_names[Dij_df$Var2], levels = rev(new_parameter_names[ord]))

# keep lower triangle incl. diagonal (each pair shown once)
Dij_df <- Dij_df[as.integer(Dij_df$Var1) + as.integer(Dij_df$Var2) <= n_par + 1, ]

# cell labels: estimate and 95% bootstrap CI; pairs whose upper CI < 0.001 shown as "<0.001"
fmt3 <- function(x) formatC(x, digits = 3, format = "f")

# larger point-estimate label
Dij_df$value_label <- ifelse(
  is.na(Dij_df$value),
  "",
  ifelse(Dij_df$hi < 0.001,
         "<0.001",
         fmt3(Dij_df$value))
)

# smaller 95% bootstrap CI label
Dij_df$ci_label <- ifelse(
  is.na(Dij_df$value) | Dij_df$hi < 0.001,
  "",
  paste0("[", fmt3(Dij_df$lo), ", ", fmt3(Dij_df$hi), "]")
)

# white text on dark cells, black text otherwise
Dij_df$txt_col <- ifelse(!is.na(Dij_df$value) & Dij_df$value > 0.6 * max(Dij, na.rm = TRUE),
                         "white", "black")

p <- ggplot(Dij_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  # point estimate, larger and slightly above the cell centre
  geom_text(
    aes(label = value_label, color = txt_col),
    size = 7,
    nudge_y = 0.11,
    fontface = "plain") +
  
  # 95% bootstrap CI, smaller and slightly below the cell centre
  geom_text(
    aes(label = ci_label, color = txt_col),
    size = 3.9,
    nudge_y = -0.13,
    fontface = "plain") +
  scale_color_identity() +
  scale_fill_gradient(low = "white", high = "deepskyblue4",
                      limits = c(0, max(Dij, na.rm = TRUE)),
                      oob = scales::squish, na.value = "grey85",
                      name = "D<sub>ij</sub>") +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.theme = element_markdown(size = 20,
                                   margin = margin(b = 15))
  )) +
  labs(x = "",
       y = "") +
  theme_light() +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1, size = 26),
        axis.text.y = element_markdown(size = 26, margin = margin(r = 8)),
        legend.title = element_markdown(size = 20),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank())

ggsave("Dij plot.svg", p, device = svglite, width = 10, height = 8, bg = "white")

# validation: analytical steady state vs. long ODE integration (L0 != R0, incl. L at the maximum)
set.seed(250724)
val_rows <- sample(seq_len(N), 200)
L_check  <- 10^c(-4, -2, 0, 2, 4)
t_grid   <- c(0, 10^seq(-4, 7, by = 0.5))   # log-spaced output times up to t = 1e7
val <- do.call(rbind, lapply(val_rows, function(i) {
  pars <- unlist(A[i, ])
  L_i  <- c(L_check, 10^max_P_over_L(pars)[["logL_opt"]])
  t(vapply(L_i, function(L0) {
    state <- c(L = L0, R = 50, C0 = 0, C1 = 0, C2 = 0, Y = 0, P = 0)
    out <- ode(y = state, times = t_grid, func = ode_system, parms = pars,
               method = "lsoda", rtol = 1e-8, atol = 1e-10, maxsteps = 1e6)
    c(P_ode   = unname(out[nrow(out), "P"]),
      P_ana   = steady_state_P(pars, L0 = L0),
      reached = as.numeric(out[nrow(out), "time"] >= max(t_grid)))
  }, numeric(3)))
}))
ok      <- val[, "reached"] == 1
rel_dev <- abs(val[ok, "P_ana"] - val[ok, "P_ode"]) / pmax(abs(val[ok, "P_ode"]), 1e-6)
cat(sprintf("Integrations reaching t = 1e7: %d of %d\n", sum(ok), length(ok)))
summary(rel_dev)
cat(sprintf("Max. relative deviation analytical vs. ODE: %.2e\n", max(rel_dev)))











#### Code part 2: Sobol screening sample size evaluation ----
# the second part calculates and plots the Sobol indix STi over different sample sizes
# afterwards it checks when stability for each parameter is achieved
# NOTE: requires steady_state_P(), evaluate_readout() and scale_to_bounds() from code part 1
# (evaluate_readout() uses the analytical steady state, no ODE integration)

rm(list = intersect(c("A", "B", "AB_list", "f_A", "f_B", "f_AB", "results"), ls()))
set.seed(250724)


#initiate results dataframe 
results <- data.frame()

time1 <- Sys.time()


# number of cores for parallel calculation
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)


results_Dij <- data.frame()
k     <- length(parameter_bounds)
pairs <- t(combn(k, 2))
pair_names <- paste(names(parameter_bounds)[pairs[, 1]], names(parameter_bounds)[pairs[, 2]], sep = " x ")


# generation of matrices for different sample sizes

for (N in c(seq(100, 1000, by = 100), seq(2000, 2^16, by = 1000), 2^16)) { # 2^16 für volle Analyse
  k <- length(parameter_bounds)
  sobol_samples <- sobol(N, 2 * k)
  A <- sobol_samples[, 1:k]
  B <- sobol_samples[, (k + 1):(2 * k)]
  
  A <- as.data.frame(scale_to_bounds(A, parameter_bounds))
  B <- as.data.frame(scale_to_bounds(B, parameter_bounds))
  colnames(A) <- colnames(B) <- names(parameter_bounds)
  
  AB_list <- lapply(1:k, function(i) {
    AB_i <- A
    AB_i[, i] <- B[, i]
    return(AB_i)
  })
  
  
  # analytical steady state for matrix A 
  f_A <- foreach(i = 1:N, .combine = c,
                 .export = c("evaluate_readout", "max_P_over_L", "steady_state_P", "logL_grid")) %dopar% {
                   evaluate_readout(A[i, ])
                 }
  
  # analytical steady state for matrix B 
  f_B <- foreach(i = 1:N, .combine = c,
                 .export = c("evaluate_readout", "max_P_over_L", "steady_state_P", "logL_grid")) %dopar% {
                   evaluate_readout(B[i, ])
                 }
  
  # analytical steady state for ABi matrices (ABi matrix = matrix A with ith column from matrix B )
  f_AB <- foreach(i = 1:k, .combine = 'c', .multicombine = TRUE, 
                  .export = c("evaluate_readout", "max_P_over_L", "steady_state_P", "logL_grid")) %dopar% {
                    list(apply(AB_list[[i]], 1, evaluate_readout))
                  }
  
  
  
  #calculate variance for the model 
  f_all <- c(f_A, f_B)
  stopifnot(all(is.finite(f_A)))
  stopifnot(all(is.finite(f_B)))
  stopifnot(all(vapply(f_AB, function(x) all(is.finite(x)), logical(1))))
  
  V_Y <- mean(f_all^2) - mean(f_all)^2

  #calculate first order Sobol index
  Si <- sapply(1:k, function(i) {
    mean(f_B * (f_AB[[i]] - f_A)) / V_Y
  })
  
  #calculate total order Sobol index
  STi <- sapply(1:k, function(i) {
    (mean((f_A - f_AB[[i]])^2) / 2) / V_Y
  })
  
  #copy results for each sample amount into temp and then save in results
  temp_result <- data.frame(N = N, Parameter = names(parameter_bounds), Si = Si, STi = STi)
  results <- rbind(results, temp_result)
  
  Dij_N <- apply(pairs, 1, function(pr) {
    STi[pr[1]] + STi[pr[2]] - (mean((f_AB[[pr[1]]] - f_AB[[pr[2]]])^2) / 2) / V_Y
  })
  results_Dij <- rbind(results_Dij, data.frame(N = N, Pair = pair_names, Dij = Dij_N))
}

stopCluster(cl)



time2 <- Sys.time()
time2-time1

# new parameter names
par_levels <- unname(new_parameter_names[names(parameter_bounds)])
results$Parameter_label <- factor(unname(new_parameter_names[results$Parameter]),
                                  levels = par_levels)

pair_labels <- paste(par_levels[pairs[, 1]], par_levels[pairs[, 2]], sep = " × ")
results_Dij$Pair_label <- factor(pair_labels[match(results_Dij$Pair, pair_names)],
                                 levels = pair_labels)



final_N <- results[results$N == max(results$N), ]
lev_Si  <- as.character(final_N$Parameter_label[order(final_N$Si,  decreasing = TRUE)])
lev_STi <- as.character(final_N$Parameter_label[order(final_N$STi, decreasing = TRUE)])


plot_Si <- ggplot(results, aes(x = N, y = Si, color = Parameter_label)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of samples",
       y = "S<sub>i</sub>",
       color = "Parameter") +
  theme_minimal() +
  theme(axis.title.y = element_markdown(),
        legend.text  = element_markdown()) +
  ylim(-0.05, NA) +
  scale_color_manual(values = pal, breaks = lev_Si)


plot_STi <- ggplot(results, aes(x = N, y = STi, color = Parameter_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = "Number of samples",
       y = "S<sub>Ti</sub>",
       color = "Parameter") +
  theme_minimal() +
  theme(
    axis.text.y  = element_markdown(size = 18, margin = margin(r = 8)),
    axis.text.x  = element_text(size = 18),
    axis.title.y = element_markdown(size = 18),
    axis.title.x = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text  = element_markdown(size = 14),
    legend.key.size = unit(1.5, "cm")
  ) +
  scale_color_manual(values = pal, breaks = lev_STi)

plot_Dij <- ggplot(results_Dij, aes(x = N, y = Dij, color = Pair_label)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Number of samples",
       y = "D<sub>ij</sub>",
       color = "Pair") +
  theme_minimal() +
  theme(axis.title.y = element_markdown(size = 18),
        legend.text  = element_markdown(size = 10))


print(plot_Dij)
print(plot_Si)
print(plot_STi)


rel_tolerance <- 0.01
abs_tolerance <- 0.005
window_size   <- 5

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

stability_Si <- do.call(
  rbind,
  lapply(unique(results$Parameter), function(parameter) {
    
    d <- results[results$Parameter == parameter, ]
    
    cbind(
      Index = "Si",
      Parameter = parameter,
      check_final_stability(
        N = d$N,
        values = d$Si,
        rel_tolerance = rel_tolerance,
        abs_tolerance = abs_tolerance,
        window_size = window_size
      )
    )
  })
)

stability_STi <- do.call(
  rbind,
  lapply(unique(results$Parameter), function(parameter) {
    
    d <- results[results$Parameter == parameter, ]
    
    cbind(
      Index = "STi",
      Parameter = parameter,
      check_final_stability(
        N = d$N,
        values = d$STi,
        rel_tolerance = rel_tolerance,
        abs_tolerance = abs_tolerance,
        window_size = window_size
      )
    )
  })
)

stability_Dij <- do.call(
  rbind,
  lapply(unique(results_Dij$Pair), function(pair) {
    
    d <- results_Dij[results_Dij$Pair == pair, ]
    
    cbind(
      Index = "Dij",
      Parameter = pair,
      check_final_stability(
        N = d$N,
        values = d$Dij,
        rel_tolerance = rel_tolerance,
        abs_tolerance = abs_tolerance,
        window_size = window_size
      )
    )
  })
)
convergence_details <- rbind(
  stability_Si,
  stability_STi,
  stability_Dij
)

all_converged <- all(convergence_details$Converged)

cat(
  "\nAll Si, STi and Dij converged:",
  all_converged,
  "\n\n"
)

print(
  convergence_details[
    ,
    c(
      "Index",
      "Parameter",
      "Converged",
      "Final_window_start_N",
      "Window_end_N",
      "Max_abs_deviation",
      "Allowed_abs_deviation"
    )
  ],
  row.names = FALSE
)
