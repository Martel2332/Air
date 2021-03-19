
alpha <- 4.48        
beta <- 0.76         
sigma2 <- 81.14      
J <- 3               
y <-c(21, 20, 15)
n <-c(48, 34, 21)
z<-c(10, 30, 50)

alpha
beta
J
y
n
z

#Les parametres a estimer sont donc x1, x2, x3, theta1 et theta2. Cependant, la loi de ej est connue, 
#Nous connaissons donc exactement la loi des x1, x2 et x3.


#Sample 10000

N=10000

x1_etoile=NULL
x2_etoile=NULL
x3_etoile=NULL


theta1=NULL
theta2=NULL
for (i in 1:N){
  
  x1_etoile[i]=rnorm(1,mean=alpha+beta*z[1],sd=sqrt(sigma2))
  
  x2_etoile[i]=rnorm(1,mean=alpha+beta*z[2],sd=sqrt(sigma2))
  
  x3_etoile[i]=rnorm(1,mean=alpha+beta*z[3],sd=sqrt(sigma2))
  
  
  
  
  
  theta1[i]= rnorm(1,0,sqrt(10000))+rnorm(1,mean=0,sd=sqrt(10000))
  
  theta2[i]=rnorm(1,0,sqrt(10000))+rnorm(1,mean=0,sd=sqrt(10000))
  
}


mean(x1_etoile)
mean(x2_etoile)
mean(x3_etoile)




my1stMHWithinGibbs <- function(nchain, mu, Sigma, init,
                               prop.sd = 0.5 * sqrt(diag(Sigma))){
  ## Just to show that an argument can be missing in a function call
  ## if you use inside the function body the missing function
  if (missing(init))
    init <- mu
  
  iSigma <- solve(Sigma)
  ## Generalization of the exercise to arbitrary dimension
  dim <- length(mu)
  
  ## Allocation of the Markov chain
  ## (please please always allocate me...
  ##  I don't want to slow down your code)
  chain <- matrix(NA, nchain + 1, dim)
  chain[1,] <- init
  
  acc.rates <- rep(0, dim)
  
  for (iter in 1:nchain){
    ## Loop over iteration
    current <- chain[iter,]
    
    for (j in 1:dim){
      ## Loop over dimension
      
      prop <- current
      prop[j] <- rnorm(1, current[j], prop.sd[j])
      
      ## Compute the acceptance proba
      top <- -0.5 * mahalanobis(prop, mu, iSigma, inverted = TRUE)
      bottom <- -0.5 * mahalanobis(current, mu, iSigma, inverted = TRUE)
      acc.prob <- exp(top - bottom)
      
      if (runif(1) < acc.prob){
        current <- prop
        acc.rates[j] <- acc.rates[j] + 1
      }
    }
    
    ## Populate the Markov chain
    chain[iter+1,] <- current      
  }
  
  return(list(chain = chain, acc.rates = acc.rates / nchain))
}

## Application

mu=0
Sigma=10000

theta1_G <- my1stMHWithinGibbs(10000, mu, Sigma)
theta2_G <- my1stMHWithinGibbs(10000, mu, Sigma)

chain1 <- theta1_G$chain
chain2=theta2_G$chain

acc.rates1 <- theta1_G$acc.rates
acc.rates1

acc.rates2=theta2_G$acc.rates
acc.rates2

par(mfrow=c(2,1))
plot(chain1)
plot(chain2)
## pair plot to check for ellipsoids typical of the multivariate normal
plot(as.data.frame(chain1))
plot(as.data.frame(chain2))


par(mfrow=c(1,5))
plot(theta2, type='l')
plot(x1_etoile, type='l')
plot(x2_etoile, type='l')
plot(x3_etoile, type='l')
plot(theta1, type='l')












