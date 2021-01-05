simulate <- function(N=1e5,
                     tmax=52*10,
                     S0=0.7,
                     I0=0.0001,
                     beta0=1.8,
                     beta1=0.15,
                     s=0.07,
                     daily.mutation=1e-4,
                     mean.mutation=0.6,
                     sd.mutation=0.4) {
  weekly.mutation <- 1-(1-daily.mutation)^7
  shape.mutation <- 1/(0.4/0.6)^2
  
  last.x <- last.y <- rep(NA, N)
  
  Ivec <- rep(0, tmax)
  Ivec[1] <- round(N*I0)
  
  statevec <- rep("S", N)
  statevec[1:round(N*I0)] <- "I"
  
  last.x[1:round(N*(1-S0))] <- rnorm(round(N*(1-S0)), 0, 0.1)
  last.y[1:round(N*(1-S0))] <- rnorm(round(N*(1-S0)), 0, 0.1)
  
  for (t in 2:tmax) {
    print(t)
    beta <- beta0 * (1 + beta1 * cos(2 * pi * t/52))
    
    winf <- which(statevec=="I")
    wsus <- which(statevec=="S")
    
    ## FOI:
    ## 1 - exp(- \sum_{i=1}^I rho_i * beta/N)
    
    for (x in wsus) {
      d.last <- sqrt((last.x[x] - last.x[winf])^2 + (last.y[x] - last.y[winf])^2)
      
      rho1 <- pmin(d.last * s, 1, na.rm=TRUE)
      
      infprob <- 1-exp(-sum(rho1 * beta)/N)
      
      if(runif(1, 0, 1) < infprob) {
        ww <- winf[sample(length(winf), 1, prob=rho1)]
        
        new.x <- last.x[ww]
        new.y <- last.y[ww]
        
        if (runif(1, 0, 1) < weekly.mutation) {
          direction <- runif(1, 0, 2*pi)
          size <- rgamma(1, shape=shape.mutation, rate=shape.mutation/mean.mutation)
          
          new.x <- new.x + size * cos(direction)
          new.y <- new.y + size * sin(direction)
          
        }
        
        last.x[x] <- new.x
        last.y[x] <- new.y
        
        statevec[x] <- "I"
      }
    }
    
    ## recovery
    statevec[winf] <- "S"
    
    Ivec[t] <- sum(statevec=="I")
  }
  
  
    
}