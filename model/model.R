
# this script contains the model function

model_function = function(t, pop, param) {
  
  with(as.list(c(param, pop)), {
    
    N=S+E1+E2+Ia+Is+R+S_c+E1_c+E2_c+Ia_c+Is_c+R_c
    
    # Community force of infection
    # estimated with a sin function, shifted by 300 to align with the start of the simulation
    lambda_v = max_lambda_v/2*(sin((2*pi/period)*t+300)+1)
    
    # Total force of infection
    # workplace infections + weekend infections + telework infections
    # For a frequency-dependent model, R0 = beta*alpha/gamma_a
    # For a density-dependent model, R0 = beta*alpha*N/gamma_a
    beta = R0*gamma_a/(prop_a)
    lambda = beta*((Ia+Ia_c)/N)*(5/7*(1-alpha)) + lambda_v*epsilon*(5/7*alpha) + lambda_v*(2/7)
    
    
    # Individuals without work-related chronic disease ####
    # Susceptible
    # - workplace infections - weekend infections - telework infections - chronic disease incidence
    dS = -lambda*S - omega*alpha*S
    
    # Exposed (non-infectious)
    # + infections - progressions to exposed infectious - chronic disease incidence
    dE1 = lambda*S - sigma*E1 - omega*alpha*E1
    
    # Exposed (infectious)
    # + progressions from exposed non-infectious - progressions to infected - chronic disease incidence
    dE2 = sigma*E1 - rho*E2 - omega*alpha*E2
    
    # Infected (asymptomatic)
    # + progressions from exposed - recoveries - chronic disease incidence
    dIa = rho*E2*prop_a - gamma_a*Ia - omega*alpha*Ia
    
    # Infected (symptomatic)
    # Note: no chronic disease incidence, as we assume these individuals are not working
    # + progressions from exposed - recoveries
    dIs = rho*E2*(1-prop_a) - gamma_s*Is
    
    # Recovered
    # + recoveries - chronic disease incidence
    dR = Ia*gamma_a + Is*gamma_s - omega*alpha*R
    
    
    
    # Individuals without work-related chronic disease ####
    # Susceptible
    # - workplace infections - weekend infections - telework infections + chronic disease incidence
    dS_c = -lambda*S_c + omega*alpha*S
    
    # Exposed (non-infectious)
    # + infections - progressions to infected + chronic disease incidence
    dE1_c = lambda*S_c - sigma*E1_c + omega*alpha*E1
    
    # Exposed (infectious)
    # + progressions from exposed non-infectious - progressions to infected + chronic disease incidence
    dE2_c = sigma*E1_c - rho*E2_c + omega*alpha*E2
    
    # Infected (asymptomatic)
    # + progressions from exposed - recoveries + chronic disease incidence
    dIa_c = rho*E2_c*prop_a - gamma_a*Ia_c + omega*alpha*Ia
    
    # Infected (symptomatic)
    # + progressions from exposed - recoveries
    dIs_c = rho*E2_c*(1-prop_a) - gamma_s*Is_c
    
    # Recovered
    # + recoveries + chronic disease incidence
    dR_c = Ia_c*gamma_a + Is_c*gamma_s + omega*alpha*R
    
    list(c(dS, dE1, dE2, dIa, dIs, dR, dS_c, dE1_c, dE2_c, dIa_c, dIs_c, dR_c))
    
  })
  
}
