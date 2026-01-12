# WORKS PERFECTLY!

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


#### Code part 1: Sobol screening for first, second and total order ----
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

# Function to find the stable value of P
find_stable_P <- function(P_values, tolerance = 1e-4, consecutive_points = 5) {
  for (i in 1:(length(P_values) - consecutive_points + 1)) {
    if (max(abs(diff(P_values[i:(i + consecutive_points - 1)]))) < tolerance) {
      return(P_values[i + consecutive_points - 1])
    }
  }
  return(NA)  # Return NA if no stable value is found
}


# number of samples (N) and parameters (k)
N <- 5000   
k <- length(parameter_bounds) 

# generate sobol sequence for matrix A and B (2 * k)
sobol_samples <- sobol(N, 2 * k, scrambling = 3)

# first k columns for A, second k columns for B
A <- sobol_samples[, 1:k]
B <- sobol_samples[, (k + 1):(2 * k)]

# Scale values to each parameter's bounds in the logarithmic space
scale_to_bounds <- function(sample, bounds) {
  mapply(function(x, range) 10^(log10(range[1]) + x * (log10(range[2]) - log10(range[1]))), 
         as.data.frame(sample), bounds)
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


# solve ODE system for each parameter set
solve_ode <- function(params) {
  state <- c(L = 50, R = 50, C0 = 0, C1 = 0, C2 = 0, Y = 0, P = 0)
  times <- seq(0, 2500, by = 1)  # simulation time
  
  out <- ode(y = state, times = times, func = ode_system, parms = params)
  
  # extract steady state value of P
  find_stable_P(out[, "P"])  
}

time1 <- Sys.time()


# number of cores for parallel calculation
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# parallel calculation for matrix A
f_A <- foreach(i = 1:N, .combine = c, .packages = 'deSolve') %dopar% {
  solve_ode(A[i, ])
}

# parallel calculation for matrix B
f_B <- foreach(i = 1:N, .combine = c, .packages = 'deSolve') %dopar% {
  solve_ode(B[i, ])
}

# parallel calculation for ABi matrices
f_AB <- foreach(i = 1:k, .combine = 'c', .multicombine = TRUE, .packages = 'deSolve') %dopar% {
  list(apply(AB_list[[i]], 1, solve_ode))
}

stopCluster(cl)

time2 <- Sys.time()
time2-time1

# calculate model variance 
f_all <- c(f_A, f_B, unlist(f_AB))
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


#second order index
Vij <- matrix(0, nrow = k, ncol = k)


# calculation according to Saltelli et al., 2010 (compare paper) 
for (i in 1:(k - 1)) {
  for (j in (i + 1):k) {
    Vij[i,j] <- (1/(2*N)) * sum((f_AB[[i]] - f_AB[[j]])^2, na.rm=TRUE) 
    Vij[j, i] <- Vij[i, j]  
  }
}

rownames(Vij) <- parameter_names
colnames(Vij) <- parameter_names

Sij <- Vij / V_Y

rownames(Sij) <- parameter_names
colnames(Sij) <- parameter_names



# total order index
STi <- sapply(1:k, function(i) {
  (mean((f_A - f_AB[[i]])^2, na.rm = TRUE) / 2) / V_Y
})

# Convert the Sij matrix to a data frame suitable for ggplot
Sij_df <- as.data.frame(as.table(Sij))

# Rename the columns to match the expected names in ggplot
colnames(Sij_df) <- c("Var1", "Var2", "value")



# new parameter names 
new_parameter_names <- c(  "kon" = "k<sub>on</sub>",
                           "koff" = "k<sub>off</sub>",
                           "kp" = "k<sub>p</sub>",
                           "phi" = "k<sub>lim</sub>",
                           "Ymax" = "M<sub>max</sub>",
                           "lambda" = "k<sub>amp</sub>",
                           "delta" = "k<sub>act</sub>",
                           "my" = "k<sub>inh</sub>")

# create df for sensitivity indices with new parameter names
result <- data.frame(Parameter = new_parameter_names[names(parameter_bounds)], Si = Si, STi = STi)

# sort df for Si
result_sorted <- result[order(result$Si, decreasing = TRUE), ]


Sij_df$Var1 <- new_parameter_names[Sij_df$Var1]
Sij_df$Var2 <- new_parameter_names[Sij_df$Var2]

# Si barplot
ggplot(result_sorted, aes(x = reorder(Parameter, Si), y = Si)) +
  geom_bar(stat = "identity", fill = "deepskyblue", color = "black") +
  labs( x = "", y = "S<sub>i</sub>") +
  theme_light() +
  coord_flip() +
  theme(
    axis.text.y = element_markdown(size = 18, margin = margin(r = 8)),     
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_markdown(size = 18),
    plot.margin = margin(10, 20, 10, 10)      
  ) +
  scale_y_continuous(
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, by = 0.2),
    expand = expansion(mult = 0, add = 0)  
  )

#sort df for STi
result_sorted_STi <- result[order(result$STi, decreasing = TRUE), ]

# STi barplot
ggplot(result_sorted_STi, aes(x = reorder(Parameter, STi), y = STi)) +
  geom_bar(stat = "identity", fill = "deepskyblue4", color = "black") +
  labs( x = "", y = "S<sub>Ti</sub>") +
  theme_light() +
  coord_flip() +
  theme(
    axis.text.y = element_markdown(size = 18, margin = margin(r = 8)),     
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_markdown(size = 18),
    plot.margin = margin(10, 20, 10, 10)       
  )+
  scale_y_continuous(
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, by = 0.2),
    expand = expansion(mult = 0, add = 0) 
  )



# create Sij (second order) indices heatmap plot and save
p <- ggplot(Sij_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "white", high = "deepskyblue4", mid = "deepskyblue", 
                       midpoint = 0.5, limit = c(0, 1), 
                       name = "Sensitivity Index") +
  guides(fill = guide_colorbar(
    title.position = "top",                          
    title.theme = element_text(size = 16, 
                               margin = margin(b = 15))  
  )) +
  labs(x = "", 
       y = "") +
  theme_light() +
  theme(axis.text.x = element_markdown(angle = 0, hjust = 1, size = 14),
        axis.text.y = element_markdown(size = 14, margin = margin(r = 8)),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),   
        legend.key.size = unit(1.5, "cm")    )    

ggsave("plot_fixed.svg", p, device = svglite, width = 8, height = 6, bg = "white")



#### Code part 2: Sobol screening sample size evaluation ----
# the second part calculates and plots the Sobol indix STi over different sample sizes
# afterwards it checks when stability for each parameter is achieved
# CAVE: in case of high sample sizes calculation will take very long!

rm(list = c("A", "B", "AB_list", "f_A", "f_B", "f_AB", "results"))  
set.seed(250724)


#initiate results dataframe 
results <- data.frame()

time1 <- Sys.time()


# number of cores for parallel calculation
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# generation of matrices for different sample sizes
for (N in c(seq(100, 1000, by = 100), seq(2000, 15000, by = 1000))) {
  k <- length(parameter_bounds)
  sobol_samples <- sobol(N, 2 * k, scrambling = 3)
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
  
  
  # solving ODEs for matrix A 
  f_A <- foreach(i = 1:N, .combine = c, .packages = 'deSolve') %dopar% {
    solve_ode(A[i, ])
  }
  
  # solving ODEs for matrix B 
  f_B <- foreach(i = 1:N, .combine = c, .packages = 'deSolve') %dopar% {
    solve_ode(B[i, ])
  }
  
  # solving ODEs for ABi matrices (ABi matrix = matrix A with ith column from matrix B )
  f_AB <- foreach(i = 1:k, .combine = 'c', .multicombine = TRUE, .packages = 'deSolve') %dopar% {
    list(apply(AB_list[[i]], 1, solve_ode))
  }
  
  
  
  #calculate variance for the model 
  f_all <- c(f_A, f_B, unlist(f_AB))
  V_Y <- mean(f_all^2, na.rm = TRUE) - mean(f_all, na.rm = TRUE)^2
  
  #calculate first order Sobol index
  Si <- sapply(1:k, function(i) {
    mean(f_B * (f_AB[[i]] - f_A), na.rm = TRUE) / V_Y
  })
  
  #calculate total order Sobol index
  STi <- sapply(1:k, function(i) {
    (mean((f_A - f_AB[[i]])^2, na.rm = TRUE) / 2) / V_Y
  })
  
  #copy results for each sample amount into temp and then save in results
  temp_result <- data.frame(N = N, Parameter = names(parameter_bounds), Si = Si, STi = STi)
  results <- rbind(results, temp_result)
}

stopCluster(cl)

time2 <- Sys.time()
time2-time1


# Plot Si over different samples sizes
plot_Si <- ggplot(results, aes(x = N, y = Si, color = Parameter)) +
  geom_line() +
  geom_point() +
  labs(x = "Anzahl der Samples", 
       y = "Si") +
  theme_minimal() +
  ylim(-0.05, NA)


# Plot STi over different sample sizes
plot_STi <- ggplot(results, aes(x = N, y = STi, color = Parameter)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Number of samples", 
    y = "S<sub>Ti</sub>") +
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(size = 18, margin = margin(r = 8)),     
    axis.text.x = element_text(size = 18),
    axis.title.y = element_markdown(size = 18),
    axis.title.x = element_text(size = 18),
    legend.title = element_text(size = 14), 
    legend.text = element_markdown(size = 14),   
    legend.key.size = unit(1.5, "cm")   
  )


print(plot_Si)
print(plot_STi)



# check for stability function
check_stability <- function(values, tolerance = 0.05) {
  stable_indices <- NULL
  for (i in 1:(length(values) - 2)) {
    relative_changes <- abs(diff(values[i:(i + 2)]) / values[i:(i + 1)])
    if (all(relative_changes < tolerance)) {
      stable_indices <- c(stable_indices, i)  # Save the first stable index
    }
  }
  return(stable_indices)
}

# initiate df
stability_summary <- data.frame(Parameter = character(), Stable_N = integer())

# check for STi stability for each parameter (5% tolerance)
for (parameter in unique(results$Parameter)) {
  param_results <- results[results$Parameter == parameter, ]
  
  # check stability
  stable_indices <- check_stability(param_results$STi, tolerance = 0.05)
  
  # when stable, save sample size (N) for the parameter
  if (length(stable_indices) > 0) {
    first_stable_N <- param_results$N[stable_indices[1]]
    stability_summary <- rbind(stability_summary,
                               data.frame(Parameter = parameter, Stable_N = first_stable_N))
  }
}

# generate text that includes stable sample size for each parameter
output_text <- paste("Stability of parameters achieved for number of samples:\n",
                     paste(stability_summary$Parameter, "-", 
                           stability_summary$Stable_N, "samples", sep=" ", collapse="\n"))

cat(output_text)
