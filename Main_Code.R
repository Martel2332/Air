Air = function(niter=10^5, init=c(0,0,0,0,0), prop.sd=c(0.4,0.015)){   #ajuster prop.sd selon acc.rates
  
  alpha = 4.48        
  beta = 0.76         
  sigma2 = 81.14
  sigma2_theta = 1024
  J = 3               
  y = c(21, 20, 15)
  n = c(48, 34, 21)
  z = c(10, 30, 50)
  
  chain=matrix(NA, niter+1, 5)
  chain[1,] = init
  
  acc.rates = rep(0, 2)
  
  for (iter in 1:niter){
    current = chain[iter,]
    
    #Mise à jour des x_i
    for (i in 1:J){
      prop = rnorm(1, alpha+beta*z[i], sqrt(sigma2))
      current[i] = prop
    }
    x = c(current[1], current[2], current[3])
    
    #Mise à jour de theta_1
    prop = rnorm(1, current[4], prop.sd[1])
    
    top = sum(y*(prop+current[5]*x) - n*log(1+exp(prop+current[5]*x)) + (-(prop^2+current[5]^2)/(2*sigma2_theta)))  
    bottom = sum(y*(current[4]+current[5]*x) - n*log(1+exp(current[4]+current[5]*x)) + (-(current[4]^2+current[5]^2)/(2*sigma2_theta)))
    acc.prob = exp(top - bottom)
    
    if (runif(1) < acc.prob){
      current[4] = prop
      acc.rates[1] = acc.rates[1] + 1
    }
    
    #Mise à jour de theta_2
    prop = rnorm(1, current[5], prop.sd[2])
    
    top = sum(y*(current[4]+prop*x) - n*log(1+exp(current[4]+prop*x)) + (-(current[4]^2+prop^2)/(2*sigma2_theta)))
    bottom = sum(y*(current[4]+current[5]*x) - n*log(1+exp(current[4]+current[5]*x)) + (-(current[4]^2+current[5]^2)/(2*sigma2_theta)))
    acc.prob = exp(top - bottom)
    
    if (runif(1) < acc.prob){
      current[5] = prop
      acc.rates[2] = acc.rates[2] + 1
    }
    
    chain[iter+1,] = current
  }
  
  return(list(chain=chain, acc.rates=acc.rates/niter))
}

acc.rates = Air()$acc.rates
chain = Air()$chain

par(mfrow=c(1,5))
for (i in 1:5){
  plot(chain[,i], type="l")
}

mean(chain[,4])
mean(chain[,5])

sd(chain[,4])
sd(chain[,5])

summary(chain[,4])
summary(chain[,5])

par(mfrow=c(1,5))
for (i in 1:5){
  acf(chain[,i], lag.max=100)
}

#A faire : Elagage des theta et peut-être burning de x3
