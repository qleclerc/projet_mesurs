
# this script runs 3 example scenarios of homeworking proportion (0, 0.5, 1)
# the community FOI is displayed in the final plot for convenience

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)

source(here::here("Model", "model.R"))

# beta = 3              # rate of infection in workplace
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

result0 = as.data.frame(lsoda(Init.cond, Time, model_function, param))
result0$alpha = 0
param["alpha"] = 0.5         # proportion of homeworking
result05 = as.data.frame(lsoda(Init.cond, Time, model_function, param), alpha = 0.5)
result05$alpha = 0.5
param["alpha"] = 1          # proportion of homeworking
result1 = as.data.frame(lsoda(Init.cond, Time, model_function, param), alpha = 1)
result1$alpha = 1

result_all = rbind(result0, result05, result1)

p1 = result_all %>%
  melt(id.vars=c("time", "alpha")) %>%
  ggplot() +
  facet_grid(cols = vars(alpha)) +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw() +
  labs(x = "Time", y = "Number of individuals", col = "")

p2 = data.frame(time = result0$time,
                val = max_lambda_v/2*(sin((2*pi/period)*result0$time+300)+1)) %>%
  ggplot() +
  geom_line(aes(time, val)) +
  theme_bw() +
  labs(x = "Time", y = "Community force of infection") 

plot_grid(p1, p2, nrow=2,
          rel_heights = c(1, 0.4))

#ggsave(here::here("figures", "example_homework_scenarios.png"))
