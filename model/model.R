
# this script contains the model function

model_function = function(t, pop, param) {
  
  with(as.list(c(param, pop)), {
    
    N=S+E+Ia+Is+R
    
    # Community force of infection
    # estimated with a sin function, shifted by 300 to align with the start of the simulation
    lambda_v = max_lambda_v/2*(sin((2*pi/period)*t+300)+1)
    
    # Total force of infection
    # workplace infections + weekend infections + homework infections
    # For a frequency-dependent model, R0 = beta*alpha/gamma_a
    # For a density-dependent model, R0 = beta*alpha*N/gamma_a
    beta = R0*gamma_a/(prop_a)
    lambda = beta*(Ia/N)*(5/7*(1-alpha)) + lambda_v*(5/7*alpha) + lambda_v*epsilon*(2/7)
    
    # Susceptible
    # - workplace infections - weekend infections - homework infections
    dS = -lambda*S
    
    # Exposed
    # + infections - progressions to infected
    dE = lambda*S - sigma*E
    
    # Infected (asymptomatic)
    # + progressions from exposed - recoveries
    dIa = sigma*E*prop_a - gamma_a*Ia
    
    # Infected (symptomatic)
    # + progressions from exposed - recoveries
    dIs = sigma*E*(1-prop_a) - gamma_s*Is
    
    # Recovered
    # + recoveries
    dR = Ia*gamma_a + Is*gamma_s
    
    list(c(dS, dE, dIa, dIs, dR))
    
  })
  
}
