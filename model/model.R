
# this script contains the model function

model_function = function(t, pop, param) {
  
  with(as.list(c(param, pop)), {
    
    N=S+E+Ia+P+Is+R+S_c+E_c+Ia_c+P_c+Is_c+R_c
    
    # Community force of infection
    # estimated with a sin function, shifted by 300 to align with the start of the simulation
    lambda_v = max_lambda_v/2*(sin((2*pi/period)*t+300)+1)
    
    # Total force of infection
    # workplace infections + weekend infections + telework infections
    # For a frequency-dependent model, R0 = beta*alpha/gamma_a
    # For a density-dependent model, R0 = beta*alpha*N/gamma_a
    
    ## TODO change beta calculation here! DONE !
    if (alpha == 1) {
      beta = 0
    } else {
      beta = R0 * rho * gamma_a / ( (1-alpha) * ((1-prop_a) * gamma_a + rho * nu * prop_a) )
    }
    lambda = 
      5/7 * (1-alpha) * beta / (N-Is) * (nu * (Ia+Ia_c) + (P+P_c)) + 
      5/7 * alpha * epsilon * lambda_v + 
      2/7 * lambda_v
    
    # Individuals without work-related chronic disease ####
    # Susceptible
    # - workplace infections - weekend infections - telework infections - chronic disease incidence
    dS = -lambda*S - omega*alpha*S
    
    # Exposed (non-infectious)
    # + infections - progressions to infectious - chronic disease incidence
    dE = lambda*S - sigma*E - omega*alpha*E
    
    # Infected (asymptomatic)
    # + progressions from exposed - recoveries - chronic disease incidence
    dIa = sigma*E*prop_a - gamma_a*Ia - omega*alpha*Ia
    
    # Pre-symptomatic (infectious)
    # + progressions from exposed - progressions to symptomatic - chronic disease incidence
    dP = sigma*E*(1-prop_a) - rho*P - omega*alpha*P
    
    # Infected (symptomatic)
    # Note: no chronic disease incidence, as we assume these individuals are not working
    # + progressions from pre-symptomatic - recoveries
    dIs = rho*P - gamma_s*Is
    
    # Recovered
    # + recoveries - chronic disease incidence
    dR = Ia*gamma_a + Is*gamma_s - omega*alpha*R
    
    
    
    # Individuals without work-related chronic disease ####
    # Susceptible
    # - workplace infections - weekend infections - telework infections + chronic disease incidence
    dS_c = -lambda*S_c + omega*alpha*S
    
    # Exposed (non-infectious)
    # + infections - progressions to infectious + chronic disease incidence
    dE_c = lambda*S_c - sigma*E_c + omega*alpha*E
    
    # Infected (asymptomatic)
    # + progressions from exposed - recoveries + chronic disease incidence
    dIa_c = sigma*E_c*prop_a - gamma_a*Ia_c + omega*alpha*Ia
    
    # Pre-symptomatic (infectious)
    # + progressions from exposed - progressions to symptomatic + chronic disease incidence
    dP_c = sigma*E_c*(1-prop_a) - rho*P_c + omega*alpha*P
    
    # Infected (symptomatic)
    # + progressions from pre-symptomatic - recoveries
    dIs_c = rho*P_c - gamma_s*Is_c
    
    # Recovered
    # + recoveries + chronic disease incidence
    dR_c = Ia_c*gamma_a + Is_c*gamma_s + omega*alpha*R
    
    list(c(dS, dE, dIa, dP, dIs, dR, dS_c, dE_c, dIa_c, dP_c, dIs_c, dR_c))
    
  })
  
}
