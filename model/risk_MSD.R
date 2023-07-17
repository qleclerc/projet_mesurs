
#daily risk 
risk_MSD = function(alpha) {
  
  
  risk = 0.18+1/exp(((exp(2.18*alpha)-1.23*exp(alpha*1.95))))
  
  return(risk/(100*52))
  
}
  