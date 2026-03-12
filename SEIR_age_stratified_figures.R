###################################################
### Clean workspace
###################################################
rm(list = ls())

#Set current working directory
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

###################################################
### Libraries
###################################################
library(deSolve)

###################################################
### Parameters
###################################################

# Total population
N0 <- 66800000

# Age-group population sizes (edit if needed)
N_B0 <- 1200000          # babies
N_C0 <- 12100000         # children
N_A0 <- N0 - N_B0 - N_C0 # adults

# Demographic rates (per day)
mu_B <- 9.7 / (1000 * 365.25)   # birth rate
mu_D <- 9.5 / (1000 * 365.25)   # death rate

# Ageing rates (per day)
a_BC <- 1 / (2 * 365.25)        # babies -> children (after 2 years on average)
a_CA <- 1 / (16 * 365.25)       # children -> adults (after 16 years on average)

# Measles transmission parameters
TIME_INFECTIOUS <- 8
TIME_LATENT <- 19
R0_MEASLES <- 18

delta <- 1 / TIME_LATENT
gamma <- 1 / TIME_INFECTIOUS
beta  <- R0_MEASLES * gamma

# Immunity quantities
adults_immune = 1-0.092
childs_immune = 0.96 #Mean value
baby_vac_rate = 0.943
vac_efficacy = 0.95

freq_vac_baby <- baby_vac_rate*vac_efficacy/(2*365.25) #in the equations  

# Pack parameters
params <- c(
  mu_B  = mu_B,
  mu_D  = mu_D,
  a_BC  = a_BC,
  a_CA  = a_CA,
  beta  = beta,
  delta = delta,
  gamma = gamma,
  freq_vac_baby = freq_vac_baby
)

###################################################
### Initial conditions
###################################################

# Initial infected / exposed
I_B0 <- 1
I_C0 <- 0
I_A0 <- 0

E_B0 <- 0
E_C0 <- 0
E_A0 <- 0

# Initial recovered
# For now, set to 0 unless you want pre-existing immunity.
R_B0 <- 0.5*N_B0 * baby_vac_rate * vac_efficacy #0.5 because vaccination after 12 months
R_C0 <- N_C0 * childs_immune
R_A0 <- N_A0 * adults_immune 

# Initial susceptibles
S_B0 <- N_B0 - E_B0 - I_B0 - R_B0
S_C0 <- N_C0 - E_C0 - I_C0 - R_C0
S_A0 <- N_A0 - E_A0 - I_A0 - R_A0

# State vector
xstart <- c(
  S_B = S_B0, E_B = E_B0, I_B = I_B0, R_B = R_B0,
  S_C = S_C0, E_C = E_C0, I_C = I_C0, R_C = R_C0,
  S_A = S_A0, E_A = E_A0, I_A = I_A0, R_A = R_A0
)

###################################################
### Simulation time
###################################################
SIMTIME <- 1000
times <- seq(0, SIMTIME, by = 0.25)

###################################################
### SEIR age-structured model
###################################################
seir_age_model <- function(t, x, params) {
  with(as.list(c(x, params)), {
    
    # Total living population at time t
    N <- S_B + E_B + I_B + R_B +
      S_C + E_C + I_C + R_C +
      S_A + E_A + I_A + R_A
    
    # Total infectious population
    Itot <- I_B + I_C + I_A
    
    # Force of infection
    lambda <- beta * Itot / N
    
    # New infections per day (instantaneous flow S -> E)
    incidence <- lambda * (S_B + S_C + S_A)
    
    ###################################################
    ### Babies
    ###################################################
    dS_B <- mu_B * N - lambda * S_B - a_BC * S_B - mu_D * S_B - freq_vac_baby * S_B
    dE_B <- lambda * S_B - a_BC * E_B - delta * E_B - mu_D * E_B
    dI_B <- delta * E_B - gamma * I_B - a_BC * I_B - mu_D * I_B
    dR_B <- gamma * I_B - a_BC * R_B - mu_D * R_B + freq_vac_baby * S_B
    
    ###################################################
    ### Children
    ###################################################
    dS_C <- -lambda * S_C + a_BC * S_B - a_CA * S_C - mu_D * S_C
    dE_C <-  lambda * S_C + a_BC * E_B - a_CA * E_C - delta * E_C - mu_D * E_C
    dI_C <-  delta * E_C + a_BC * I_B - a_CA * I_C - gamma * I_C - mu_D * I_C
    dR_C <-  gamma * I_C + a_BC * R_B - a_CA * R_C - mu_D * R_C
    
    ###################################################
    ### Adults
    ###################################################
    dS_A <- -lambda * S_A + a_CA * S_C - mu_D * S_A
    dE_A <-  lambda * S_A + a_CA * E_C - delta * E_A - mu_D * E_A
    dI_A <-  delta * E_A + a_CA * I_C - gamma * I_A - mu_D * I_A
    dR_A <-  gamma * I_A + a_CA * R_C - mu_D * R_A
    
    list(c(
      dS_B, dE_B, dI_B, dR_B,
      dS_C, dE_C, dI_C, dR_C,
      dS_A, dE_A, dI_A, dR_A
    ),
    incidence = incidence)
  })
}

###################################################
### Helper functions
###################################################

# Create initial state for a given scenario
make_initial_state <- function(baby_vac_rate, childs_immune, vac_efficacy,
                               N0, N_B0, N_C0, adults_immune,
                               E_B0, E_C0, E_A0, I_B0, I_C0, I_A0) {
  
  N_A0 <- N0 - N_B0 - N_C0
  
  R_B0 <- 0.5 * N_B0 * baby_vac_rate * vac_efficacy
  R_C0 <- N_C0 * childs_immune
  R_A0 <- N_A0 * adults_immune
  
  S_B0 <- N_B0 - E_B0 - I_B0 - R_B0
  S_C0 <- N_C0 - E_C0 - I_C0 - R_C0
  S_A0 <- N_A0 - E_A0 - I_A0 - R_A0
  
  c(
    S_B = S_B0, E_B = E_B0, I_B = I_B0, R_B = R_B0,
    S_C = S_C0, E_C = E_C0, I_C = I_C0, R_C = R_C0,
    S_A = S_A0, E_A = E_A0, I_A = I_A0, R_A = R_A0
  )
}

# Create parameter vector for a given scenario
make_params <- function(mu_B, mu_D, a_BC, a_CA, beta, delta, gamma,
                        baby_vac_rate, vac_efficacy) {
  
  freq_vac_baby <- baby_vac_rate * vac_efficacy / (2 * 365.25)
  
  c(
    mu_B  = mu_B,
    mu_D  = mu_D,
    a_BC  = a_BC,
    a_CA  = a_CA,
    beta  = beta,
    delta = delta,
    gamma = gamma,
    freq_vac_baby = freq_vac_baby
  )
}

# Run one scenario and return daily output
run_scenario <- function(baby_vac_rate, childs_immune) {
  
  params_run <- make_params(
    mu_B = mu_B, mu_D = mu_D, a_BC = a_BC, a_CA = a_CA,
    beta = beta, delta = delta, gamma = gamma,
    baby_vac_rate = baby_vac_rate, vac_efficacy = vac_efficacy
  )
  
  xstart_run <- make_initial_state(
    baby_vac_rate = baby_vac_rate,
    childs_immune = childs_immune,
    vac_efficacy = vac_efficacy,
    N0 = N0, N_B0 = N_B0, N_C0 = N_C0,
    adults_immune = adults_immune,
    E_B0 = E_B0, E_C0 = E_C0, E_A0 = E_A0,
    I_B0 = I_B0, I_C0 = I_C0, I_A0 = I_A0
  )
  
  out <- as.data.frame(
    lsoda(y = xstart_run, times = times, func = seir_age_model, parms = params_run)
  )
  
  dt <- diff(times)[1]
  out$newInfections_step <- out$incidence * dt
  out$day <- floor(out$time)
  
  daily_incidence <- aggregate(newInfections_step ~ day, data = out, sum)
  out.daily <- out[out$time %% 1 == 0, ]
  out.daily$newInfections <- daily_incidence$newInfections_step[
    match(out.daily$time, daily_incidence$day)
  ]
  
  out.daily$I_total <- out.daily$I_B + out.daily$I_C + out.daily$I_A
  out.daily$E_total <- out.daily$E_B + out.daily$E_C + out.daily$E_A
  out.daily$S_total <- out.daily$S_B + out.daily$S_C + out.daily$S_A
  out.daily$R_total <- out.daily$R_B + out.daily$R_C + out.daily$R_A
  
  out.daily
}

# Save all plots for one scenario
save_scenario_plots <- function(out.daily, folder, baby_vac_rate, childs_immune) {
  
  baby_vac_label  <- round(baby_vac_rate, 2)
  child_imm_label <- round(childs_immune, 2)
  
  title_suffix <- paste0(
    "(baby_vac = ", baby_vac_label,
    ", child_immun = ", child_imm_label, ")"
  )
  
  file_suffix <- paste0(
    "_babyVac_", formatC(baby_vac_rate, format = "f", digits = 3),
    "_childImm_", formatC(childs_immune, format = "f", digits = 3),
    ".png"
  )
  
  ymax_seir <- max(
    out.daily$S_total,
    out.daily$E_total,
    out.daily$I_total,
    out.daily$R_total,
    na.rm = TRUE
  )
  
  ymax_I_age <- max(
    out.daily$I_B, out.daily$I_C, out.daily$I_A,
    na.rm = TRUE
  )
  
  # 1) New infected
  png(file.path(folder, paste0("new_infected", file_suffix)), width = 1400, height = 900, res = 150)
  plot(
    out.daily$time, out.daily$newInfections,
    type = "l", lwd = 3, col = 2,
    xlab = "Time (days)", ylab = "New infections per day",
    main = paste("New Infected", title_suffix)
  )
  dev.off()
  
  # 2) SEIR
  png(file.path(folder, paste0("SEIR", file_suffix)), width = 1400, height = 900, res = 150)
  plot(
    out.daily$time, out.daily$S_total,
    type = "l", lwd = 3, col = 2,
    ylim = c(0, ymax_seir),
    xlab = "Time (days)", ylab = "Population",
    main = paste("SEIR", title_suffix)
  )
  lines(out.daily$time, out.daily$E_total, lwd = 3, col = 6)
  lines(out.daily$time, out.daily$I_total, lwd = 3, col = 3)
  lines(out.daily$time, out.daily$R_total, lwd = 3, col = 4)
  legend(
    "right",
    legend = c("S", "E", "I", "R"),
    col = c(2, 6, 3, 4),
    lwd = 3,
    bty = "n"
  )
  dev.off()
  
  # 3) Infected by group
  png(file.path(folder, paste0("infected_by_group", file_suffix)), width = 1400, height = 900, res = 150)
  plot(
    out.daily$time, out.daily$I_B,
    type = "l", lwd = 3, col = 2,
    ylim = c(0, ymax_I_age),
    xlab = "Time (days)", ylab = "Infectious",
    main = paste("Infected by Group", title_suffix)
  )
  lines(out.daily$time, out.daily$I_C, lwd = 3, col = 3)
  lines(out.daily$time, out.daily$I_A, lwd = 3, col = 4)
  legend(
    "right",
    legend = c("I_B", "I_C", "I_A"),
    col = c(2, 3, 4),
    lwd = 3,
    bty = "n"
  )
  dev.off()
  
  # 4) Total infected
  png(file.path(folder, paste0("total_infected", file_suffix)), width = 1400, height = 900, res = 150)
  plot(
    out.daily$time, out.daily$I_total,
    type = "l", lwd = 3, col = 4,
    xlab = "Time (days)", ylab = "Total infectious",
    main = paste("Total Infected", title_suffix)
  )
  dev.off()
}

# Save summary plots for a list of scenarios
save_summary_plots <- function(results_list, folder, varying_name) {
  
  scenario_names <- names(results_list)
  
  # Summary: new infected
  ymax_new <- max(sapply(results_list, function(x) max(x$newInfections, na.rm = TRUE)))
  png(file.path(folder, paste0("summary_new_infected_", varying_name, ".png")),
      width = 1400, height = 900, res = 150)
  
  plot(
    results_list[[1]]$time, results_list[[1]]$newInfections,
    type = "l", lwd = 2,
    ylim = c(0, ymax_new),
    xlab = "Time (days)", ylab = "New infections per day",
    main = paste("Summary New Infected -", varying_name)
  )
  
  if (length(results_list) > 1) {
    for (i in 2:length(results_list)) {
      lines(results_list[[i]]$time, results_list[[i]]$newInfections, lwd = 2, col = i)
    }
  }
  
  legend(
    "topright",
    legend = scenario_names,
    col = seq_along(results_list),
    lwd = 2,
    bty = "n",
    cex = 0.9
  )
  dev.off()
  
  # Summary: total infected
  ymax_total <- max(sapply(results_list, function(x) max(x$I_total, na.rm = TRUE)))
  png(file.path(folder, paste0("summary_total_infected_", varying_name, ".png")),
      width = 1400, height = 900, res = 150)
  
  plot(
    results_list[[1]]$time, results_list[[1]]$I_total,
    type = "l", lwd = 2,
    ylim = c(0, ymax_total),
    xlab = "Time (days)", ylab = "Total infectious",
    main = paste("Summary Total Infected -", varying_name)
  )
  
  if (length(results_list) > 1) {
    for (i in 2:length(results_list)) {
      lines(results_list[[i]]$time, results_list[[i]]$I_total, lwd = 2, col = i)
    }
  }
  
  legend(
    "topright",
    legend = scenario_names,
    col = seq_along(results_list),
    lwd = 2,
    bty = "n",
    cex = 0.9
  )
  dev.off()
}

###################################################
### Create figure folders
###################################################

figures_dir <- file.path(getwd(), "Figures")
dir.create(figures_dir, showWarnings = FALSE)

folder_baby_vac <- file.path(figures_dir, "Figures_baby_vac_rate")
folder_child_immun <- file.path(figures_dir, "Figures_childs_immune")
folder_both <- file.path(figures_dir, "Figures_baby_vac_rate_and_childs_immune")

dir.create(folder_baby_vac, showWarnings = FALSE)
dir.create(folder_child_immun, showWarnings = FALSE)
dir.create(folder_both, showWarnings = FALSE)

###################################################
### Parameter grids
###################################################

baby_vac_values <- c(0.0, 0.3, 0.6, 0.9, 0.95, 0.975, 0.99, 1.0)
child_immun_values <- c(0.0, 0.3, 0.6, 0.9, 0.96, 0.975, 0.99, 1.0)
both_values <- c(0.0, 0.3, 0.6, 0.9, 0.95, 0.975, 0.99, 1.0)

###################################################
### 1) Vary baby_vac_rate only
###################################################

results_baby_vac <- list()

for (val in baby_vac_values) {
  out.daily <- run_scenario(
    baby_vac_rate = val,
    childs_immune = childs_immune
  )
  
  results_baby_vac[[paste0("baby_vac=", formatC(val, format = "f", digits = 3))]] <- out.daily
  
  save_scenario_plots(
    out.daily = out.daily,
    folder = folder_baby_vac,
    baby_vac_rate = val,
    childs_immune = childs_immune
  )
}

save_summary_plots(
  results_list = results_baby_vac,
  folder = folder_baby_vac,
  varying_name = "baby_vac_rate"
)

###################################################
### 2) Vary childs_immune only
###################################################

results_child_immun <- list()

for (val in child_immun_values) {
  out.daily <- run_scenario(
    baby_vac_rate = baby_vac_rate,
    childs_immune = val
  )
  
  results_child_immun[[paste0("child_immun=", formatC(val, format = "f", digits = 3))]] <- out.daily
  
  save_scenario_plots(
    out.daily = out.daily,
    folder = folder_child_immun,
    baby_vac_rate = baby_vac_rate,
    childs_immune = val
  )
}

save_summary_plots(
  results_list = results_child_immun,
  folder = folder_child_immun,
  varying_name = "childs_immune"
)

###################################################
### 3) Vary both together
###################################################

results_both <- list()

for (val in both_values) {
  out.daily <- run_scenario(
    baby_vac_rate = val,
    childs_immune = val
  )
  
  results_both[[paste0("baby_vac=child_immun=", formatC(val, format = "f", digits = 3))]] <- out.daily
  
  save_scenario_plots(
    out.daily = out.daily,
    folder = folder_both,
    baby_vac_rate = val,
    childs_immune = val
  )
}

save_summary_plots(
  results_list = results_both,
  folder = folder_both,
  varying_name = "baby_vac_rate_and_childs_immune"
)