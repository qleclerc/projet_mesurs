
# this script allows you to quickly run the model and plot the output

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)

source(here::here("Model", "model.R"))

#beta = 3             # rate of infection in workplace
R0 = 3.3              # Basic reproduction number
alpha = 0             # proportion of homeworking
max_lambda_v = 0.01   # maximum community force of infection
period = 200          # duration of epidemic wave in community
epsilon = 0.5         # relative force of infection during homeworking
sigma = 1/3           # progression rate from exposed to infected
prop_a = 0.2          # proportion of asymptomatic infections
gamma_a = 1/5         # recovery rate for asymptomatics
gamma_s = 1/5         # recovery rate for symptomatics

dt = 0.1    # time-step
Tmax = 200  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,Is=0,R=0) 
param = c(beta, max_lambda_v, period, epsilon, prop_a, gamma_a, gamma_s, sigma)

result = as.data.frame(lsoda(Init.cond, Time, model_function, param))

result %>%
  melt(id.vars="time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw() +
  labs(x = "Time", y = "Number of individuals", col = "")

#ggsave(here::here("figures", "example_simple_run.png"))