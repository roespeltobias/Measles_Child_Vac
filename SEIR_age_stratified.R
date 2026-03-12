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
childs_immune = 1#0.96 #Mean value
baby_vac_rate = 1#0.943
vac_efficacy = 1#0.95

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
    ))
  })
}

###################################################
### Solve the system
###################################################
out <- as.data.frame(lsoda(y = xstart, times = times, func = seir_age_model, parms = params))

###################################################
### Daily output
###################################################
out.daily <- out[out$time %% 1 == 0, ]

# Total infectious / exposed / susceptible / recovered
out.daily$I_total <- out.daily$I_B + out.daily$I_C + out.daily$I_A
out.daily$E_total <- out.daily$E_B + out.daily$E_C + out.daily$E_A
out.daily$S_total <- out.daily$S_B + out.daily$S_C + out.daily$S_A
out.daily$R_total <- out.daily$R_B + out.daily$R_C + out.daily$R_A

# Approximate daily new infections from decline in S_total
out.daily$newInfections <- c(
  out.daily$I_total[1],
  out.daily$S_total[-nrow(out.daily)] - out.daily$S_total[-1]
)

###################################################
### Inspect results
###################################################
head(out.daily)
tail(out.daily)

max(out.daily$newInfections)
out.daily$time[which.max(out.daily$newInfections)]

###################################################
### Plots
###################################################

# Daily new infections
plot(
  out.daily$time, out.daily$newInfections,
  type = "l", lwd = 3, col = 2,
  xlab = "Time (days)", ylab = "New infections per day"
)

# Total SEIR
ymax_seir <- max(
  out.daily$S_total,
  out.daily$E_total,
  out.daily$I_total,
  out.daily$R_total,
  na.rm = TRUE
)

plot(
  out.daily$time, out.daily$S_total,
  type = "l", lwd = 3, col = 2,
  ylim = c(0, ymax_seir),
  xlab = "Time (days)", ylab = "Population"
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

# Infectious by age group
ymax_I <- max(out.daily$I_B, out.daily$I_C, out.daily$I_A, na.rm = TRUE)

plot(
  out.daily$time, out.daily$I_B,
  type = "l", lwd = 3, col = 2,
  ylim = c(0, ymax_I),
  xlab = "Time (days)", ylab = "Infectious"
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