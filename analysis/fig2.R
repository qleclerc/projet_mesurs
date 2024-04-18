
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

R0 = 2.66              # Basic reproduction number
alpha = 0.1           # proportion of teleworking
t_alpha = -1          # activation time for teleworking (default -1 = always on)
nu = 0.35             # relative force of infection of asymptomatic cases
epsilon = 0.5         # relative force of infection during teleworking
sigma = 1/6.57         # progression rate from exposed to infectious
rho = 1/1.5           # progression rate from pre-symptomatic to symptomatic
prop_a = 0.2          # proportion of asymptomatic infections
gamma_a = 1/5         # recovery rate for asymptomatics
gamma_s = 1/5         # recovery rate for symptomatics
baseline_NCD = 0.0005   # baseline NCD rate

dt = 0.1    # time-step
Tmax = 90  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,P=0,Is=0,R=0,S_c=0,E_c=0,Ia_c=0,P_c=0,Is_c=0,R_c=0) 
param = c(R0=R0, alpha=alpha, t_alpha=t_alpha,
          nu=nu, epsilon=epsilon, sigma=sigma, rho=rho, prop_a=prop_a,
          gamma_a=gamma_a, gamma_s=gamma_s, baseline_NCD=baseline_NCD)

# community force of infection
# uses approxfun() to generate an interpolating function, passed to the model function
# whilst there is no data, still using sin function shifted by 300 to align with start of simulation
commu_FOI = approxfun(c(0:1000), 0.005/2*(sin((2*pi/90)*c(0:1000)+300)+1))

#with data will look something like:
# commu_FOI = approxfun(data_time_points, data_values)

# rate of chronic disease linked to telework
# uses approxfun() to generate an interpolating function, passed to the model function
# whilst there is no data, still using a constant max rate multiplied by alpha

result_NC_all = data.frame()
result_RR_NCD = data.frame(alpha = alpha_to_try)

# LINEAR INCREASING ####

NCD_rate = approxfun(c(0,0.5,1), c(7.04/7.04,7.3/7.04,7.46/7.04))
DRF = "LI"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}

result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))

result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# LINEAR DECREASING ####

# NCD_rate = approxfun(c(0,1/20,4/20,1), c(-3.33/-3.33, -0.571/3.33+1, -0.997/3.33+1, -1.233/3.33+1))
NCD_rate = approxfun(c(0,1), c(-3.33/-3.33, -1.233/3.33+1))
DRF = "LD"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}

result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))
result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# INVERTED U SHAPED ####

NCD_rate = approxfun(c(0,0.2,0.5,1), c(49.9/49.9, 52.8/49.9, 53.2/49.9, 47.5/49.9))
DRF = "IU"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}

result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))
result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# U SHAPED ####

NCD_rate = approxfun(c(0,0.5,1), c(2.33/2.33, 1.71/2.33, 2.08/2.33))
DRF = "US"

result_all = data.frame()

for(alpha_t in alpha_to_try){
  param["alpha"] = alpha_t         # proportion of teleworking
  result_all = rbind(result_all,
                     data.frame(lsoda(Init.cond, Time, model_function, param,
                                      commu_FOI = commu_FOI, NCD_rate = NCD_rate),
                                alpha = alpha_t))
}

result_NC_all = rbind(result_NC_all,
                      result_all %>% filter(time == max(time)) %>% mutate(DRF = DRF))

result_NC_all = result_NC_all %>%
  mutate(DRF = factor(DRF, levels = c("LI", "LD", "US", "IU")))
result_RR_NCD[,DRF] = NCD_rate(alpha_to_try)

# PLOT ####

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
  geom_line(aes(time, value, colour = variable), linewidth = 1) +
  scale_colour_brewer(palette = "Dark2") +
  theme_bw() +
  labs(x = "Time (days)", y = "Number of individuals", col = "") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow=1))

p2 = result_NC_all %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c) %>%
  ggplot() +
  geom_point(aes(DRF, Tot_c, shape = DRF), size = 3) +
  facet_grid(cols = vars(alpha)) +
  theme_bw() +
  guides(shape = "none") +
  labs(x = "Dose-response curve", y = "Cumulative number of individuals\neventually developing a NCD", col = "")


plot_grid(p1, p2,
          nrow=2,
          labels = c("a)", "b)"), hjust = 0)

ggsave(here::here("figures", "fig2.png"), height = 7, width = 7)

melt(result_RR_NCD, id.vars = c("alpha")) %>%
  rename(DRF=variable, RR=value) %>%
  mutate(DRF=factor(DRF, levels = c("LI", "LD", "US", "IU"))) %>%
  ggplot() +
  geom_line(aes(alpha, RR, colour = DRF), linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  labs(x = "Telework proportion", y = "Relative risk of developing a NCD", col = "") +
  scale_y_continuous(limits = c(0.6, 1.1)) +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=11))

ggsave(here::here("figures", "figpres.png"))

