# remove all object from the working space 
rm(list = ls())
# set working directory
setwd('C:\\Users\\sikde\\OneDrive\\Studium\\M.Sc. Statistics\\5_semester\\Analysis')


#function of tree simulation
simTree_anti <- function(b, d, S0, sigma, T, r) {
  tot <- sum(b^(1:(d - 1)))
  S1 <- numeric(tot + 1)
  S2 <- numeric(tot + 1)
  S1[1] <- S0
  S2[1] <- S0
  dt <- T/d
  for (i in 0:(tot - b^(d - 1))) {
    for (j in 1:b) {
      z =  rnorm(1)
      #stock price for z
      S1[i * b + j + 1] <- S1[i + 1] *
        exp((r - 0.5 * sigma^2) *
              dt + sigma * sqrt(dt) *z)
      #stock price for -z
      S2[i * b + j + 1] <- S2[i + 1] *
        exp((r - 0.5 * sigma^2) *
              dt + sigma * sqrt(dt) *(-z))
    }
  }
  return(list(S1, S2))
}



#function of upper estimator
upperBG <- function(S, b, d, f, r) {
  tot <- sum(b^(1:(d - 1)))
  start <- tot - b^(d - 1) + 1
  end <- tot + 1
  P<-S
  P[start:end] <- f(S[start:end])
  tot1 <- sum(b^(1:(d - 2)))
  for (i in tot1:0) {
    m <- mean(P[i * b + 1:b + 1])*exp(-r)
    v <- f(S[i + 1])
    P[i + 1] <- max(v, m)
  }
  P
}



#function lower estimator
lowerBG <- function(S, b, d, f, r) {
  tot <- sum(b^(1:(d - 1)))
  start <- tot - b^(d - 1) + 1
  end <- tot + 1
  p<-S
  p[start:end] <- f(S[start:end])
  tot1 <- sum(b^(1:(d - 2)))
  m <- numeric(b)
  for (i in tot1:0) {
    v <- f(S[i + 1])
    for (j in 1:b) {
      m[j] <- mean(p[i * b + (1:b)[-j] + 1])*exp(-r)
      m[j] <- ifelse(v > m[j], v, p[i * b + (1:b)[j] +
                                      + 1]*exp(-r))
    }
    p[i + 1] <- mean(m)
  }
  p
}


#payoff function
f <- function(x) sapply(x, function(x) max(x - K, 0))


#function of Confidence interval with respect to upper and lower estimator
CI = function(f, low_est, upper_est, alpha, n, S0){
  lower_interval = max(f(S0), (mean(low_est) - qnorm(1-alpha/2)*(sd(low_est)/sqrt(n))))
  higher_interval = mean(upper_est) + qnorm(1-alpha/2)*(sd(upper_est)/sqrt(n))
  
  interval = c(lower_interval, higher_interval)
  names(interval) = c('lower', 'upper')
  return(interval)
}


Total_time = list()
for(j in 1:7){
start = Sys.time()
#initialize Parameter value 
b <- 3; d <- j; T <- 1; r <- 0.05; sigma <- 0.4
S0 <- 36; K = 30; n = 10000; alpha = 0.05



#variable  initialisation
low_est1 = vector('numeric', n)
upper_est1 = vector('numeric', n)
low_est2 = vector('numeric', n)
upper_est2 = vector('numeric', n)
C1 = vector('numeric', n)
C2 = vector('numeric', n)


for(i in 1:n){
  #tree simulation
  S <- simTree_anti(b, d, S0, sigma, T, r)
  
  #lower estimator for i-th replication
  low_est1[i] = lowerBG(S[[1]], b, d, f, r)[1]
  low_est2[i] = lowerBG(S[[2]], b, d, f, r)[1]
  #upper estimator for i-th replication
  upper_est1[i] = upperBG(S[[1]], b, d, f, r)[1]
  upper_est2[i] = upperBG(S[[2]], b, d, f, r)[1]
  
  #weighted estimator
  C1[i] = 0.5*max(f(S0),low_est1[i]) + 0.5*upper_est1[i]
  C2[i] = 0.5*max(f(S0),low_est2[i]) + 0.5*upper_est2[i]
  
}
end = Sys.time()
Total_time[[j]] = end - start
}


# estimators with antithetic variates
low_est_anti = (low_est1 + low_est2)/2
upper_est_anti = (upper_est1 + upper_est2)/2
C_anti = (C1 + C2)/2

#Calculate the standard deviation of the estimators
#lower estimator
sd(low_est1)
sd(low_est_anti)

#upper estimator
sd(upper_est1)
sd(upper_est_anti)



#simple simulation MC estimator
MC_est_low = mean(low_est1)
MC_est_upper = mean(upper_est1)
MC_C = mean(C1)

#antithetic MC estimator
MC_low_est_anti = mean(low_est_anti)
MC_upper_est_anti = mean(upper_est_anti)
MC_C_anti = mean(C_anti)


#compute simple simulation CI interval with respect to upper and lower estimator
CI_est = CI(f, low_est1, upper_est1, alpha, n, S0)
#compute AV CI interval  with respect to upper and lower estimator
CI_anti =  CI(f, low_est_anti, upper_est_anti, alpha, n, S0)

#CI of the weighted estimator for simple MC simulation
MC_C_lower = mean(C1) - qnorm(1-alpha/2)*(sd(C1)/sqrt(n))
MC_C_upper = mean(C1) + qnorm(1-alpha/2)*(sd(C1)/sqrt(n))

#CI of the weighted estimator for MC simulation with AV
MC_C_lower = mean(C_anti) - qnorm(1-alpha/2)*(sd(C_anti)/sqrt(n))
MC_C_upper = mean(C_anti) + qnorm(1-alpha/2)*(sd(C_anti)/sqrt(n))

#result simple simulation
result = c(MC_est_low, MC_est_upper,  MC_C, sd(low_est1)/sqrt(n), 
           sd(upper_est1)/sqrt(n), sd(C1)/sqrt(n), CI_est)
names(result) = c('Lower est.', 'Upper est.', 'C', 'sigma lower', 
                  'sigma upper', 'sigma C', 'CI lower', 'CI upper')
result_est = round(result, digits = 2)
#result_est[c(3,6,7,8)]


#result antithetic simulation
result = c(MC_low_est_anti, MC_upper_est_anti, MC_C_anti, 
           sd(low_est_anti)/sqrt(n), sd(upper_est_anti)/sqrt(n), sd(C_anti)/sqrt(n), CI_anti)
names(result) = c('Lower est.', 'Upper est.', 'C', 'sigma lower', 
                  'sigma upper', 'sigma C', 'CI lower', 'CI upper')
result_anti = round(result, digits = 2)
#result_anti[c(3,6,7,8)]

data = t(as.data.frame(result_est))
res_final = rbind(data, result_anti)
row.names(res_final) = NULL
res_final = round(res_final, digits = 2)

#export result as tex data
library(stargazer)
stargazer(res_final, type = 'latex', digits = 2)




#computation cost analysis
final_time = as.numeric(Total_time)
for(i in 1:4){
  final_time[i] = final_time[i]/60
}
final_time1 = round(final_time, digits = 2)
a1 = c(0.007998943, 0.016003132, 0.007999897, 
       0.008002043, 0.008010864, 0.007996082, 0.000000000)
a2 = c(0.015993834, 0.016013145, 0.007999897, 
       0.008002996, 0.008002996, 0.016003132, 0.000000000)
final_time2 = a1 + a2
final_time2 = final_time2/60
time = c(final_time1, final_time2)
method = c(rep('Random Tree', 7), rep('LSM', 7))
method = as.factor(method)
x_axis = c(1:7, 1:7)
x_axis = as.factor(x_axis)

data_time = data.frame(method, x_axis, time)
pdf(file = "C:\\Users\\sikde\\OneDrive\\Studium\\M.Sc. Statistics\\5_semester\\time.pdf",width = 11, height = 5)
colors = c('lightcoral', 'indianred4', 'lightskyblue4',
            'midnightblue', 'olivedrab', 'mediumaquamarine',
           'violetred4')

ggplot(data_time, aes(x=method, y=time, fill=x_axis)) + 
  geom_bar(stat="identity", position="dodge") +
  theme(legend.text = element_text(size = 30),
        text = element_text(size=30),
        legend.key.size = unit(0.80, 'cm'))+
  ylab("Execution time") +
  guides(fill=guide_legend(title="Period (d)"))+
  scale_fill_manual(values=colors) + theme(panel.background = element_blank())
dev.off()
