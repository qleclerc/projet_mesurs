library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

SIR.model <- function(t, pop, param) {
  with(as.list(c(param, pop)), {
    
    N=S+E+Ia+Is+R
    
    dS = -beta*(Ia/N)*S - lambda_v*S - lambda_v*epsilon*S
    dE = beta*(Ia/N)*S + lambda_v*S + lambda_v*epsilon*S - sigma*E
    dIa = sigma*E*prop_a - gamma_a*Ia
    dIs = sigma*E*(1-prop_a) - gamma_s*Is
    dR = Ia*gamma_a + Is*gamma_s
    
    list(c(dS, dE, dIa, dIs, dR))
  })
}

beta = 5
lambda_v = 0.1
epsilon = 0.5
prop_a = 0.2
gamma_a = 1/5
gamma_s = 1/5
sigma = 1/3

dt = 0.1
Tmax = 100
N = 5000

I0 = 1
S0 = N-I0

Time = seq(from=0,to=Tmax,by=dt)
Init.cond = c(S=S0,E=0,Ia=I0,Is=0,R=0) 
param = c(beta, lambda_v, epsilon, prop_a, gamma_a, gamma_s, sigma)

result = as.data.frame(lsoda(Init.cond, Time, SIR.model, param))

result %>%
  melt(id.vars="time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw()
