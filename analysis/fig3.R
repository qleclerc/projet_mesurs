
# the community FOI is displayed in the final plot for convenience

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(ggh4x)

source(here::here("Model", "model.R"))

R0 = 3.3              # Basic reproduction number
alpha = 0.8             # proportion of teleworking
t_alpha = -1          # activation time for teleworking (default -1 = always on)
nu = 0.35             # relative force of infection of asymptomatic cases
epsilon = 0.5         # relative force of infection during teleworking
sigma = 1/1.5         # progression rate from exposed to infectious
rho = 1/1.5           # progression rate from pre-symptomatic to symptomatic
prop_a = 0.2          # proportion of asymptomatic infections
gamma_a = 1/5         # recovery rate for asymptomatics
gamma_s = 1/5         # recovery rate for symptomatics

dt = 0.1    # time-step
Tmax = 90  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,P=0,Is=0,R=0,S_c=0,E_c=0,Ia_c=0,P_c=0,Is_c=0,R_c=0) 
param = c(R0=R0, alpha=alpha,
          nu=nu, epsilon=epsilon, sigma=sigma, rho=rho, prop_a=prop_a,
          gamma_a=gamma_a, gamma_s=gamma_s)

# community force of infection
commu_FOI = approxfun(c(0:1000), 0.001/2*(sin((2*pi/90)*c(0:1000)+300)+1))

result_all = data.frame()


for(alpha in seq(0.2, 1, 0.2)){
  
  param["alpha"] = alpha
  
  for(t_alpha in c(0:90)){
    
    param["t_alpha"] = t_alpha
    
    # LINEAR INCREASING ####
    
    chronic_rate = approxfun(seq(0,1,0.1), (0.01+seq(0,1,0.1)*0.06)/365)
    DRF = "LI"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, chronic_rate = chronic_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # LINEAR DECREASING ####

    chronic_rate = approxfun(seq(0,1,0.1), (0.07-seq(0,1,0.1)*0.06)/365)
    DRF = "LD"

    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, chronic_rate = chronic_rate))

    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # U SHAPED ####
    
    chronic_rate = approxfun(seq(0,1,0.1), 0.00015/2*(sin((2*pi/1)*seq(0,1,0.1)+1.5)+1.1))
    DRF = "US"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, chronic_rate = chronic_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # INVERTED U SHAPED ####

    chronic_rate = approxfun(seq(0,1,0.1), 0.00015/2*(sin((2*pi/1)*seq(0,1,0.1)+4.7)+1.1))
    DRF = "IU"

    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, chronic_rate = chronic_rate))

    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
  }
  
}

# PLOT ####

result_all2 = result_all %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
         Tot_r = R + R_c) %>%
  group_by(DRF) %>%
  mutate(Tot_c = Tot_c/max(Tot_c)*100,
         Tot_r = Tot_r/max(Tot_r)*100) %>%
  ungroup %>%
  select(DRF, alpha, t_alpha, Tot_c, Tot_r) %>% 
  mutate(Tot_tot = Tot_c+Tot_r) %>%
  mutate(DRF = factor(DRF, levels = c("LI", "LD", "US", "IU")))

low_thresholds = result_all2 %>%
  filter(Tot_c <= 50 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == min(t_alpha))
high_thresholds = result_all2 %>%
  filter(Tot_c <= 50 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == max(t_alpha))
thresholds = result_all2 %>%
  filter(Tot_c <= 50 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == max(t_alpha) | t_alpha == min(t_alpha)) %>%
  ungroup %>%
  group_by(DRF, alpha) %>%
  summarise(min_t = min(t_alpha), max_t = max(t_alpha))

ggplot(result_all2) +
  geom_rect(data = thresholds, aes(xmin = min_t, xmax = max_t, ymin=0, ymax=100), alpha = 0.2) +
  geom_line(aes(t_alpha, Tot_c, colour = "Non-communicable disease"), linewidth = 0.8, alpha = 0.3) +
  geom_line(aes(t_alpha, Tot_r, colour = "Infectious disease"), linewidth = 0.8, alpha = 0.3) +
  geom_line(data = result_all2 %>% filter(Tot_c <= 50), aes(t_alpha, Tot_c, colour = "Non-communicable disease"), linewidth = 1) +
  geom_line(data = result_all2 %>% filter(Tot_r <= 50), aes(t_alpha, Tot_r, colour = "Infectious disease"), linewidth = 1) +
  geom_point(data = result_all2 %>% filter(Tot_c == 100), aes(t_alpha, Tot_c), size = 2) +
  facet_nested_wrap(~DRF+alpha, ncol = 5) +
  scale_x_continuous(breaks = seq(0,90,15)) +
  theme_bw() +
  labs(x = "Day of teleworking implementation start", y = "Relative cumulative incidence (%)", colour = "") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11))

ggsave(here::here("figures", "fig3.png"), width = 9)

       