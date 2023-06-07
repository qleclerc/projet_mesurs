
# this script runs example scenarios of teleworking proportion
# you can define the scenarios to try by changing this vector:
alpha_to_try = c(0, 0.2, 0.4, 0.6, 0.8, 1)

# the community FOI is displayed in the final plot for convenience

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)

source(here::here("Model", "model.R"))

R0 = 3.3              # Basic reproduction number
alpha = 0             # proportion of teleworking
omega = 1/100/365     # maximum rate of chronic disease due to teleworking
max_lambda_v = 0.0001   # maximum community force of infection
period = 200          # duration of epidemic wave in community
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
param = c(R0=R0, alpha=alpha, omega=omega, max_lambda_v=max_lambda_v, period=period,
          nu=nu, epsilon=epsilon, sigma=sigma, rho=rho, prop_a=prop_a,
          gamma_a=gamma_a, gamma_s=gamma_s)

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all, data.frame(lsoda(Init.cond, Time, model_function, param), alpha = alpha_t))
}

p1 = result_all %>%
  mutate(
    S = S + S_c, 
    E = E + E_c, 
    Ia = Ia + Ia_c,
    P = P + P_c,
    Is = Is + Is_c,
    R = R + R_c
  ) %>%
  select(!contains("_c")) %>%
  melt(id.vars=c("time", "alpha")) %>%
  ggplot() +
  facet_grid(cols = vars(alpha)) +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw() +
  labs(x = "Time (days)", y = "Number of individuals", col = "")

p2 = result_all %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c) %>%
  ggplot() +
  geom_line(aes(time, Tot_c)) +
  facet_grid(cols = vars(alpha)) +
  theme_bw() +
  labs(x = "Time (days)", y = "Cumulative number of individuals eventually\ndeveloping a chronic disease", col = "")

p3 = data.frame(time = result_all$time[result_all$alpha == 0],
                val = max_lambda_v/2*(sin((2*pi/period)*result_all$time[result_all$alpha == 0]+300)+1)) %>%
  ggplot() +
  geom_line(aes(time, val)) +
  theme_bw() +
  labs(x = "Time", y = "Community force of infection") 

plot_grid(p1, p2, p3, nrow=3,
          rel_heights = c(1, 1, 0.5))

#ggsave(here::here("figures", "example_telework_scenarios.png"))
