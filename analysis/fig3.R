
# the community FOI is displayed in the final plot for convenience

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(ggh4x)

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
baseline_NCD = 0.005   # baseline NCD rate

dt = 0.1    # time-step
Tmax = 90  # max time
N = 5000    # population size

I0 = 0      # initial workplace infected
S0 = N-I0   # initial workplace susceptibles

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,P=0,Is=0,R=0,S_c=0,E_c=0,Ia_c=0,P_c=0,Is_c=0,R_c=0) 
param = c(R0=R0, alpha=alpha,
          nu=nu, epsilon=epsilon, sigma=sigma, rho=rho, prop_a=prop_a,
          gamma_a=gamma_a, gamma_s=gamma_s, baseline_NCD=baseline_NCD)

# community force of infection
commu_FOI = approxfun(c(0:1000), 0.005/2*(sin((2*pi/90)*c(0:1000)+300)+1))

result_all = data.frame()


for(alpha in seq(0, 1, 0.2)){
  
  param["alpha"] = alpha
  
  for(t_alpha in c(0:90)){
    
    param["t_alpha"] = t_alpha
    
    # LINEAR INCREASING ####
    
    # NCD_rate = approxfun(seq(0,1,0.1), (0.01+seq(0,1,0.1)*0.06)/365)
    NCD_rate = approxfun(c(0,0.5,1), c(7.04/7.04,7.3/7.04,7.46/7.04))
    DRF = "LI"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # LINEAR DECREASING ####

    # NCD_rate = approxfun(seq(0,1,0.1), (0.07-seq(0,1,0.1)*0.06)/365)
    NCD_rate = approxfun(c(0,1/20,4/20,1), c(-3.33/-3.33, -0.571/3.33+1, -0.997/3.33+1, -1.233/3.33+1))
    DRF = "LD"

    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))

    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # U SHAPED ####
    
    # NCD_rate = approxfun(seq(0,1,0.1), 0.00015/2*(sin((2*pi/1)*seq(0,1,0.1)+1.5)+1.1))
    NCD_rate = approxfun(c(0,0.5,1), c(2.33/2.33, 1.71/2.33, 2.08/2.33))
    DRF = "US"
    
    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))
    
    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
    # INVERTED U SHAPED ####

    # NCD_rate = approxfun(seq(0,1,0.1), 0.00015/2*(sin((2*pi/1)*seq(0,1,0.1)+4.7)+1.1))
    NCD_rate = approxfun(c(0,0.2,0.5,1), c(49.9/49.9, 52.8/49.9, 53.2/49.9, 47.5/49.9))
    DRF = "IU"

    res = data.frame(lsoda(Init.cond, Time, model_function, param,
                           commu_FOI = commu_FOI, NCD_rate = NCD_rate))

    result_all = rbind(result_all,
                       res %>% filter(time == max(time)) %>% mutate(DRF = DRF, alpha = alpha, t_alpha = t_alpha))
    
  }
  
}

# PLOT ####

result_all2 = result_all %>%
  mutate(Tot_c = S_c + E_c + Ia_c + P_c + Is_c + R_c,
         Tot_r = R + R_c) %>%
  group_by(DRF) %>%
  mutate(Tot_c = Tot_c/first(Tot_c)*100,
         Tot_r = Tot_r/first(Tot_r)*100) %>%
  ungroup %>%
  select(DRF, alpha, t_alpha, Tot_c, Tot_r) %>% 
  mutate(Tot_tot = Tot_c+Tot_r) %>%
  mutate(DRF = factor(DRF, levels = c("LI", "LD", "US", "IU")))


low_thresholds = result_all2 %>%
  filter(Tot_c <= 90 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == min(t_alpha))
high_thresholds = result_all2 %>%
  filter(Tot_c <= 90 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == max(t_alpha))
thresholds = result_all2 %>%
  filter(Tot_c <= 90 & Tot_r <= 50) %>%
  group_by(DRF, alpha) %>%
  filter(t_alpha == max(t_alpha) | t_alpha == min(t_alpha)) %>%
  ungroup %>%
  group_by(DRF, alpha) %>%
  summarise(min_t = min(t_alpha), max_t = max(t_alpha))

pa = ggplot(result_all2 %>% filter(alpha!=0)) +
  geom_tile(aes(x = t_alpha, y = alpha, fill = Tot_c)) +
  scale_fill_distiller(palette = "RdBu") +
  facet_wrap(~DRF) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  scale_x_continuous(breaks = seq(0,90,15)) +
  labs(x = "Day of teleworking implementation start", y = "Telework frequency", fill = "ID relative\ncumulative incidence")
  
pb = ggplot(result_all2 %>% filter(alpha!=0) %>% filter(DRF == "LI")) +
  geom_tile(aes(x = t_alpha, y = alpha, fill = Tot_r)) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  scale_x_continuous(breaks = seq(0,90,15)) +
  labs(x = "Day of teleworking implementation start", y = "Telework frequency", fill = "NCD relative\ncumulative incidence")

plot_grid(pa, pb, ncol = 1)


ggplot(result_all2) +
  geom_rect(data = thresholds, aes(xmin = min_t, xmax = max_t, ymin=0, ymax=110), alpha = 0.2) +
  geom_line(aes(t_alpha, Tot_c, colour = "Non-communicable disease"), linewidth = 0.8, alpha = 0.3) +
  geom_line(aes(t_alpha, Tot_r, colour = "Infectious disease"), linewidth = 0.8, alpha = 0.3) +
  geom_line(data = result_all2 %>% filter(Tot_c <= 90), aes(t_alpha, Tot_c, colour = "Non-communicable disease"), linewidth = 1) +
  geom_line(data = result_all2 %>% filter(Tot_r <= 50), aes(t_alpha, Tot_r, colour = "Infectious disease"), linewidth = 1) +
  facet_nested_wrap(~DRF+alpha, ncol = 6) +
  scale_x_continuous(breaks = seq(0,90,15)) +
  scale_y_continuous(breaks = seq(10,110,20), limits = c(0,110)) +
  theme_bw() +
  labs(x = "Day of teleworking implementation start", y = "Relative cumulative incidence (%)", colour = "") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11))

ggsave(here::here("figures", "fig3.png"), width = 9, height = 8)

       