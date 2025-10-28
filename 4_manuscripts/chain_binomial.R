
time_points <- 10
a <- 1
beta <- 0.25
y <- 0.5

S1 <- vector(length = time_points + 1)
S2 <- vector(length = time_points + 1)

I1 <- vector(length = time_points + 1)
I2 <- vector(length = time_points + 1)

R1 <- vector(length = time_points + 1)
R2 <- vector(length = time_points + 1)

S1[1] <- 100 
S2[1] <- 50

I1[1] <- 1
I2[1] <- 1


for (t in 2:(time_points + 1)) {
  I1[t] <- rbinom(n = 1, size = S1[t-1], p = (I1[t-1] + I2[t-1] * a)*beta)
  I2[t] <- rbinom(n = 1, size = S2[t-1], p = (I1[t-1] + I2[t-1] * a)*beta*y)
  
  S1[t] <- S1[t-1] - I1[t]
  S2[t] <- S2[t-1] - I2[t]
  
  
}