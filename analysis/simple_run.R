
# this script allows you to quickly run the model and plot the output

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)

source(here::here("Model", "model.R"))

R0 = 3.3              # Basic reproduction number
alpha = 0             # proportion of teleworking
omega = 0.001         # maximum rate of chronic disease due to teleworking
max_lambda_v = 0.01   # maximum community force of infection
period = 200          # duration of epidemic wave in community
epsilon = 0.5         # relative force of infection during teleworking
sigma = 1/1.5         # progression rate from exposed to exposed infectious
rho = 1/1.5           # progression rate from exposed infectious to infected
prop_a = 0.2          # proportion of asymptomatic infections
gamma_a = 1/5         # recovery rate for asymptomatics
gamma_s = 1/5         # recovery rate for symptomatics

dt = 0.1    # time-step
Tmax = 200  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E1=0,E2=0,Ia=I0,Is=0,R=0,S_c=0,E1_c=0,E2_c=0,Ia_c=0,Is_c=0,R_c=0) 
param = c(R0, alpha, omega, max_lambda_v, period, epsilon, sigma, rho, prop_a, gamma_a, gamma_s)

result = as.data.frame(lsoda(Init.cond, Time, model_function, param))

result %>%
  melt(id.vars="time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw() +
  labs(x = "Time", y = "Number of individuals", col = "")

#ggsave(here::here("figures", "example_simple_run.png"))