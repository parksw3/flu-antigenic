odefun <- function(t, y, param) {
  with(param, {
    
    if (istart <= t & t <= iend) {
      reduction <- (1-npi)
    } else {
      reduction <- 1
    }
    
    beta <- R0 * (1 + phi * cos(2*pi*t)) * gamma * reduction
    
    S <- y[1:n]
    I <- y[(n+1):(2*n)]
    
    deriv <- c(
      mu - beta * S * (sigmamat %*% I) - mu * S,
      beta * S * I - gamma * I - m * c(head(I, -1), 0) + m * c(0, head(I, -1)) - mu * I
    )
    
    list(deriv)
  })
}

simulate_ode <- function(R0=1.8,
                         phi=0.15,
                         gamma=52,
                         mu=1/80,
                         m=0.1,
                         d=7,
                         n=300,
                         tmax=25,
                         istart=20,
                         iend=20.6,
                         npi=0.2) {
  tvec <- seq(0, tmax, by=1/300)
  
  sigmamat <- outer(1:n, 1:n, function(x, y) exp(-((x-y)/d)^2))
  
  parms <- list(
    R0=R0,
    phi=phi,
    gamma=gamma, 
    mu=mu,
    m=m,
    n=n,
    sigmamat=sigmamat,
    istart=istart,
    iend=iend,
    npi=npi
  )
  
  y0 <- c(1-1e-6 ,rep(1, n-1), 1e-6, rep(0, n-1))
  names(y0)[1:n] <- paste0("S", 1:n)
  names(y0)[(n+1):(2*n)] <- paste0("I", 1:n)
  
  out <- as.data.frame(deSolve::ode(y0, tvec, odefun, parms))
  
  out
}
