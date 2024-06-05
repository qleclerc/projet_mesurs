
#daily risk 

NCD_function = function(DRF) {
  
  if(DRF == "LS"){
    
  NCD_rate = approxfun(c(0, 0.25/5, 1/5, 5/5),
                       c((1/(1+exp(-(-3.33)))) / (1/(1+exp(-(-3.33)))),
                         (1/(1+exp(-(-3.33-0.571)))) / (1/(1+exp(-(-3.33)))),
                         (1/(1+exp(-(-3.33-0.997)))) / (1/(1+exp(-(-3.33)))),
                         (1/(1+exp(-(-3.33-1.233)))) / (1/(1+exp(-(-3.33))))))

  } else if(DRF == "US"){
    NCD_rate = approxfun(c(0,0.5,1), c(2.33/2.33, 1.71/2.33, 2.08/2.33))
  } else if(DRF == "IU"){
    NCD_rate = approxfun(c(0,0.2,0.5,1), c(49.9/49.9, 52.8/49.9, 53.2/49.9, 47.5/49.9))
  } else stop("Supported DRF are: LS, US, IU")
  
  return(NCD_rate)

}
  