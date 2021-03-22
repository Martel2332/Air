Air = function(niter=10^5, init=c(12,27,43,-0.5,0.03), prop.sd=c(0.4,0.015)){   #ajuster prop.sd selon acc.rates
  
  alpha = 4.48        
  beta = 0.76         
  sigma2 = 81.14
  sigma2_theta = 1000
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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Etude
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

acc.rates = Air()$acc.rates
chain = Air()$chain

par(mfrow=c(1,5))
para = c("X[1]","X[2]","X[3]","theta[1]","theta[2]")
for (i in 1:5){
  plot(chain[,i], type="l", xlab="Iterations", ylab=para[i])
}

mean(chain[,1])
mean(chain[,2])
mean(chain[,3])
mean(chain[,4])
mean(chain[,5])

sd(chain[,4])
sd(chain[,5])

summary(chain[,4])
summary(chain[,5])

par(mfrow=c(1,5))
for (i in 1:5){
  acf(chain[,i], lag.max=100, ylab=para[i])
}

#Trop de corrélations pour les theta -----> élagage 
chain_elag = chain[seq(1, nrow(chain), by = 50),]

par(mfrow=c(2,5))
for (i in 1:5){
  plot(chain_elag[,i], type="l", xlab="Iterations", ylab=para[i])
  acf(chain_elag[,i], lag.max=1000, ylab=para[i])
}

mean(chain_elag[,1])
mean(chain_elag[,2])
mean(chain_elag[,3])
mean(chain_elag[,4])
mean(chain_elag[,5])

#Confrontation numérique vs théorie

alpha = 4.48        
beta = 0.76         
sigma2 = 81.14
sigma2_theta = 1000
J = 3               
y = c(21, 20, 15)
n = c(48, 34, 21)
z = c(10, 30, 50)

for (i in 1:3){
  hist(chain_elag[,i], freq=FALSE, main="", xlab=para[i])
  lines(density(rnorm(-20:80, alpha+beta*z[i], sqrt(sigma2))), lwd=2, col="blue")
}





for (i in 4:5){
  hist(chain_elag[,i], freq=FALSE, main="", xlab=para[i])
  lines(exp(sum(y*(theta1+theta2*rnorm(-20:80, alpha+beta*z[i], sqrt(sigma2))) 
               - n*log(1+exp(theta1+theta2*rnorm(-20:80, alpha+beta*z[i], sqrt(sigma2)))) 
               + (-(theta1^2+theta2^2)/(2*sigma2_theta)))), lwd=2, col="blue")
}



theta1 = theta2 = rnorm(100, 0, sqrt(sigma2_theta))

exp(sum(y*(theta1+theta2*rnorm(100, alpha+beta*z[1], sqrt(sigma2))) 
        - n*log(1+exp(theta1+theta2*rnorm(-20:80, alpha+beta*z[1], sqrt(sigma2)))) 
        + (-(theta1^2+theta2^2)/(2*sigma2_theta))))
