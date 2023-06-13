
# this script allows you to quickly run the model and plot the output

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(ggpubr)

source(here::here("Model", "model.R"))

R0 = 3.3              # Basic reproduction number
alpha = 0.1           # proportion of teleworking
omega = 1/100/365     # maximum rate of chronic disease due to teleworking
nu = 0.35             # relative force of infection of asymptomatic cases
epsilon = 0.5         # relative force of infection during teleworking
sigma = 1/1.5         # progression rate from exposed to infectious
rho = 1/1.5           # progression rate from pre-symptomatic to symptomatic
prop_a = 0.2          # proportion of asymptomatic infections
gamma_a = 1/5         # recovery rate for asymptomatics
gamma_s = 1/5         # recovery rate for symptomatics


dt = 0.1    # time-step
Tmax = 200  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,P=0,Is=0,R=0,S_c=0,E_c=0,Ia_c=0,P_c=0,Is_c=0,R_c=0) 
param = c(R0=R0, alpha=alpha, omega=omega,
          nu=nu, epsilon=epsilon, sigma=sigma, rho=rho, prop_a=prop_a,
          gamma_a=gamma_a, gamma_s=gamma_s)

# community force of infection
# uses approxfun() to generate an interpolating function, passed to the model function
# whilst there is no data, still using sin function shifted by 300 to align with start of simulation
commu_FOI = approxfun(c(0:1000), 0.001/2*(sin((2*pi/200)*c(0:1000)+300)+1))

#with data will look something like:
# commu_FOI = approxfun(data_time_points, data_values)


result = as.data.frame(lsoda(Init.cond, Time, model_function, param, commu_FOI = commu_FOI))

p1 = result %>%
  mutate(
    S = S + S_c, 
    E = E + E_c, 
    Ia = Ia + Ia_c,
    P = P + P_c,
    Is = Is + Is_c,
    R = R + R_c
    ) %>%
  select(!contains("_c")) %>%
  melt(id.vars="time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time (days)", y = "", col = "", title = "Number of individuals")

p2 = result %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c) %>%
  ggplot() +
  geom_line(aes(time, Tot_c)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time (days)", y = "", col = "", title = "Number of individuals that will\ndevelop a chronic disease")

ggarrange(p1, p2, common.legend = T, legend = "bottom", align = "hv")

#ggsave(here::here("figures", "example_simple_run.png"))