remove(list=ls())
set.seed(5)

S0 <- 36
K <- 30
T <- 1
r <- 0.05
sigma <- 0.4
d<-7
n<-10000

RG <- function (n){ 
  #Generates 2 random normally distributed vectors for every path
  # in order to generate AV. 
  z1<-rnorm(n)
  z2<- -z1
  # Both sequences are stored in the same vector. 
  z<-c(z1,z2)
  return(z)
}

LSM <- function(n, d, S0, K, sigma, r, T) { 
  
  #Vector that saves operation time
  
  Total_time <<- vector('numeric', 7)
  start <- Sys.time()
  
  #Calculates ST
  
  
  s0 <- S0/K  #Transformation (1 if S0=K)
  dt <- T/d   #Calculate delta t
  
  z <- RG(n) #Generate random normal variable for each path
  z_paths <<- data.frame(data=z )  #Save the z-path for every t. t: 1...T
  
  #Calculate value of underlying at time T 
  # Brownian Motion - Black Scholes Model
  
  # AV (1...n)
  s.t <- s0 * exp((r - 1/2 * sigma^2) * T + sigma * z * (T^0.5)) 
  # AV (n+1...2n)
  s.t[(n + 1):(2 * n)] <- s0 * exp((r - 1/2 * sigma^2) * T -
                                     sigma * z * (T^0.5)) 
  
  #Calculates Payoff Function at Time T
  
  CC <- pmax(s.t-1, 0) 
 
  
  #Returns payoff function in present value
  
  # Generates AV for payoff function <-
  # (payoff1 + payoff2)/2
  
  payoffeu <- exp(-r * T) * (CC[1:n] + CC[(n + 1):(2 * n)])/2 *K 
  
  
 
  
  #Payoff estimation at time T (European Option)
  
  euprice <- mean(payoffeu) 
  
  end = Sys.time()
  Total_time[d] <<- end - start
 
  for (k in (d - 1):1) { 
    
    start <- Sys.time()
    
    #Regressive function from sT-1 to s1
    
    z <- RG(n) #Generate random normal variable for each path
    z_paths <<- data.frame(z, z_paths) #Save generated path 
    
    #Calculate st using a Brownian bridge in order
    # to reduce computational effort
    
    #AV (1...n)
    
    #Calculate Drift between st and st+1
    mean <- (log(s0) + k * log(s.t[1:n]))/(k + 1) 
    # Calculate Volatility  between St and St+1
    vol <- (k * dt/(k + 1))^0.5 * z
    #Calculates st
    s.t.1 <- exp(mean + sigma * vol) 
    
    # AV (n+1...2n)
    # Same procedure as above
    mean <- (log(s0) + k * log(s.t[(n + 1):(2 * n)]))/(k + 1) 
    s.t.1[(n + 1):(2 * n)] <- exp(mean - sigma * vol)
    
    #Calculate payoff function a time t
    CE <- pmax(s.t.1-1, 0) 
    #Identify paths that are in the money at time t
    #This paths are selected for the regression
    #Store the indices
    
    idx <- (1:(2 * n))[CE > 0] 
    
    #Future cash flows are discounted using 
    # the free rate interest 
    
    discountedCC <- CC[idx] * exp(-r * dt) 
    
    #The basis functions for the regression are defined as constant 
    #and the first three  Laguerre polynomials
    
    basis1 <- exp(-s.t.1[idx]/2)
    basis2 <- basis1 * (1 - s.t.1[idx])
    basis3 <- basis1 * (1 - 2 * s.t.1[idx] + (s.t.1[idx]^2)/2)
    
    #Generate linear regression model
    
    #Dependent Variable - Discounted payoff >0
    #Independent Variables - Laguerre Polynomials
    
    p <- lm(discountedCC ~ basis1 + basis2 + basis3)$coefficients 
    
    #Estimate Continuity Values
    
    estimatedCC <- p[1] + p[2] * basis1 + p[3] * basis2 +
      p[4] * basis3
    
    #Comparison between estimated continuation values and the payoff function
    
    EF <- rep(0, 2 * n)
    EF[idx] <- (CE[idx] > estimatedCC)
    
    # The continuation value is composed by two parts:
    #  EF = 1  continuation value of early exercise of the option
    #  (Exercise the option on that time)
    #  EF = 0 continuation value for future time step
    #  (Wait to exercise the option)
    
    # Calculate the continuation value
    
    CC <- (EF == 0) * CC * exp(-r * dt) + (EF == 1) * CE
    s.t <- s.t.1
    
    end = Sys.time()
    Total_time[k] <<- end - start
  }
  
  start <- Sys.time()
  
  #Calculate payoff for each path
  # AV <- (payoff1+payoff2/2)
  
  payoff <- exp(-r * dt) * (CC[1:n] + CC[(n + 1):(2 * n)])/2
  
  # Estimate payoff of the option
  
  usprice <- mean(payoff * K)
  
  #Calculate error
  std <- sd(payoff * K)/sqrt(n)
  error <- 1.96 * std
  U_bound <- usprice + error
  L_bound <- usprice - error
  
  #Calculate difference between us price with eu price
  
  difference <- usprice - euprice
  colnames(z_paths) <<- seq(1,d,1)
  print(data.frame(usprice,std,U_bound, L_bound))
  
  total_time <- end - start 
  Total_time <<- Total_time + total_time
  
}


LSM(10000, d, S0, K, sigma, r, T)
 


