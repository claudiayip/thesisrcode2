aE <- 10
bE <- 2
alphad <- matrix(1,1,7)
# dirichlet for time of the day
alphan <- matrix(1,1,48)

aZ0 <- 200
bZ0 <- 200 
cZ0 <- 200
aZ1 <- 200
bZ1 <- 200
cZ1 <- 200
aZ2 <- 200
bZ2 <- 200
cZ2 <- 200


MMPP3 <- function(N, aL, bL,sig){
  aE <- as.numeric(quantile(N, probs = 0.9))
  bE <- 2
  # dirichlet for day of the week
  alphad <- matrix(sum(N)/7,1,7)
  # dirichlet for time of the day
  alphan <- matrix(sum(N)/336,7,48)
  lambda0 <- mean(N) #rgamma(1,aL, bL)
  set.seed(204)
  delta <- rdirichlet(1,alphad)*7
  set.seed(205)
  eta <- matrix(0,7,48)
  for (j in 1:7){
    eta[j,] <- rdirichlet(1, alphan[j,] )*48
  }
  lambdat <- lambda0 * eta * as.vector(delta)
  er <-  array(rep(0,48*7*5), dim = c(48,7,5))
  er[,,c(1:5)] <- t(lambdat)
  lambdatv <- as.vector(er)
  
  A0 <- matrix(1/3,3,3)
  aLS <- aL ; bLT <- bL
  abcdZ1 <- matrix(1680/9,1,3); abcdZ2 <- matrix(1680/9,1,3)
  abcdZ3 <- matrix(1680/9,1,3); abcdZ4 <- matrix(1680/9,1,3)
  alphadsj <- alphad
  alphansij <- alphan
  lambda0v <- c()
  zv <- cbind()
  for (kk in 1:150) {
    poisbinl <- function(tt){
      b <- as.numeric(0)
      c <- as.numeric(0)
      for (i in 0:N[tt]){
        b[i+1] <- dpois((N[tt]-i), lambdatv[tt])
        c[i+1] <- dnbinom(i, aE, bE/(1+bE))
      }
      return(b%*%c)
    }
    
    
    Pntzt <- function(tt,z){
      if (N[tt]==0){
        if (z == 0){likelihood <- 1}
        if (z == 1){likelihood <- exp(-lambdatv[tt])}
        if (z == 2){likelihood <- 0}
      } else {
        if (z == 0){likelihood <- 0}
        if (z == 1){likelihood <- dpois(N[tt], lambdatv[tt])}
        if (z == 2){likelihood <- poisbinl(tt)}
      }
      return(likelihood)
    }
    
    
    pznznm10 <- A0[,1]
    pznznm11 <- A0[,2]
    pznznm12 <- A0[,3]
    alphazn1v<- matrix(0,3,1680)
    alphazn1v[,1] <- c(1/3,1/3,1/3)
    for (t in 2:1680){
      if (N[t]==0){
        a1 <- sum(alphazn1v[,(t-1)]* pznznm10)
        a2 <- exp(-lambdatv[t])*sum(alphazn1v[,(t-1)]* pznznm11)
        a3 <- 0
        alphazn1v[1,t] <- a1/(a1+a2+a3)
        alphazn1v[2,t] <- a2/(a1+a2+a3)
        alphazn1v[3,t] <- a3/(a1+a2+a3)
      } else {
        a1 <- 0
        a2 <- dpois(N[t], lambdatv[t]) * sum(alphazn1v[,(t-1)]* pznznm11)
        a3 <- poisbinl(t) * sum(alphazn1v[,(t-1)]* pznznm12)
        alphazn1v[1,t] <- a1/(a1+a2+a3)
        alphazn1v[2,t] <- a2/(a1+a2+a3)
        alphazn1v[3,t] <- a3/(a1+a2+a3)
      }
    }
    
    
    
    
    
    betazN <- 1/3
    betazn1v <- matrix(0,3,1680)
    betazn1v[,1680] <- betazN
    for (t in 1679:1){
      if (N[t+1]==0) {
        b1 <- betazn1v[1,t+1] * A0[1,1] +
          betazn1v[2,t+1] * exp(-lambdatv[t+1]) * A0[1,2]
        b2 <- betazn1v[1,t+1] * A0[2,1] + 
          betazn1v[2,t+1] * exp(-lambdatv[t+1]) * A0[2,2]
        b3 <- betazn1v[1,t+1] * A0[2,1] + 
          betazn1v[2,t+1] * exp(-lambdatv[t+1]) * A0[3,2]
        betazn1v[1,t] <- b1/(b1+b2+b3)
        betazn1v[2,t] <- b2/(b1+b2+b3)
        betazn1v[3,t] <- b3/(b1+b2+b3)
      } else {
        b1 <- betazn1v[2,t+1] * dpois(N[t+1], lambdatv[t+1]) * A0[1,2] +
          betazn1v[3,t+1]* poisbinl(t+1) * A0[1,3]
        b2 <- betazn1v[2,t+1] * dpois(N[t+1], lambdatv[t+1]) * A0[2,2] +
          betazn1v[3,t+1]* poisbinl(t+1) * A0[2,3]
        b3 <- betazn1v[2,t+1] * dpois(N[t+1], lambdatv[t+1]) * A0[3,2] +
          betazn1v[3,t+1]* poisbinl(t+1) * A0[3,3]
        betazn1v[1,t] <- b1/(b1+b2+b3)
        betazn1v[2,t] <- b2/(b1+b2+b3)
        betazn1v[3,t] <- b3/(b1+b2+b3)
        
      }
    }
    
    gamma <- matrix(0,3,1680)
    g1 <- alphazn1v[1,]*betazn1v[1,] 
    g2 <- alphazn1v[2,]*betazn1v[2,] 
    g3 <- alphazn1v[3,]*betazn1v[3,]
    gamma[1,] <- g1 / (g1+g2+g3)
    gamma[2,] <- g2 / (g1+g2+g3)
    gamma[3,] <- g3 / (g1+g2+g3)
    
    pztzt1 <- matrix(0,9,1679)
    for (t in 1:1679){
      pzt0zt10 <- alphazn1v[1,t]  * Pntzt(t+1,0) * A0[1,1] *betazn1v[1,t+1]
      pzt1zt10 <- alphazn1v[2,t]  * Pntzt(t+1,0) * A0[2,1] *betazn1v[1,t+1]
      pzt2zt10 <- alphazn1v[3,t]  * Pntzt(t+1,0) * A0[3,1] *betazn1v[1,t+1]
      pzt0zt11 <- alphazn1v[1,t]  * Pntzt(t+1,1) * A0[1,2] *betazn1v[2,t+1]
      pzt1zt11 <- alphazn1v[2,t]  * Pntzt(t+1,1) * A0[2,2] *betazn1v[2,t+1]
      pzt2zt11 <- alphazn1v[3,t]  * Pntzt(t+1,1) * A0[3,2] *betazn1v[2,t+1]
      pzt0zt12 <- alphazn1v[1,t]  * Pntzt(t+1,2) * A0[1,3] *betazn1v[3,t+1]
      pzt1zt12 <- alphazn1v[2,t]  * Pntzt(t+1,2) * A0[2,3] *betazn1v[3,t+1]
      pzt2zt12 <- alphazn1v[3,t]  * Pntzt(t+1,2) * A0[3,3] *betazn1v[3,t+1]
      
      p1 <- pzt0zt10 + pzt1zt10 + pzt2zt10 # p(zt1 = 0)
      p2 <- pzt0zt11 + pzt1zt11 + pzt2zt11 # p(zt1 = 1)
      p3 <- pzt0zt12 + pzt1zt12 + pzt2zt12 # p(zt1 = 2)
      if (p1 == 0) {p1 = 1}
      if (p2 == 0) {p2 = 1}
      if (p3 == 0) {p3 = 1}
      pztzt1[1,t] <- pzt0zt10 / p1 # p(zt = 0 | zt1 = 0)
      pztzt1[2,t] <- pzt1zt10 / p1 # p(zt = 1 | zt1 = 0)
      pztzt1[3,t] <- pzt2zt10 / p1 # p(zt = 2 | zt1 = 0)
      pztzt1[4,t] <- pzt0zt11 / p2 # p(zt = 0 | zt1 = 1)
      pztzt1[5,t] <- pzt1zt11 / p2 # p(zt = 1 | zt1 = 1)
      pztzt1[6,t] <- pzt2zt11 / p2 # p(zt = 2 | zt1 = 1)
      pztzt1[7,t] <- pzt0zt12 / p3 # p(zt = 0 | zt1 = 2)
      pztzt1[8,t] <- pzt1zt12 / p3 # p(zt = 1 | zt1 = 2)
      pztzt1[9,t] <- pzt2zt12 / p3 # p(zt = 2 | zt1 = 2)
    }
    
    #print(pztzt1[,c(1:5)])
    zt <- as.numeric(0)
    zt[1680] <- sample(c(0,1,2),1, prob = gamma[,1680])
    for (tt in 1679:1){
      if (zt[tt+1] ==0 ){zt[tt] <- sample(c(0,1,2), 1,prob = pztzt1[c(1,2,3),tt])}
      else if (zt[tt+1] ==1 ){zt[tt] <- sample(c(0,1,2), 1,prob = pztzt1[c(4,5,6),tt])}
      else {zt[tt] <- sample(c(0,1,2), 1,prob = pztzt1[c(7,8,9),tt])}
    }
    print(table(zt))
    
    
    samplePNB <- function(tt){
      if (N[tt]>0){
        poisnbin <- function(i){
          return(dpois(N[tt] - i, lambdatv[tt])* dnbinom(i,aE, bE/(1+bE)))}
        pbv <- as.numeric(0)
        for (i in 0:(N[tt])){
          pbv[i+1] <- poisnbin(i)
        }
        relpbv <- pbv/sum(pbv)
        No <- sample(c(0:N[tt]),1, prob= relpbv)
      }
      else {No <- 0}
      return(No)
    }
    
    NO <- as.numeric(0)
    Ne <- as.numeric(0)
    for (j in 1:1680){
      if (zt[j]==0){
        NO[j] <- 0
        Ne[j] <- 0
      }
      if (zt[j]==1){
        NO[j] <- N[j]
        Ne[j] <- 0
      } else if (zt[j] == 2){
        NO[j] <- samplePNB(j)
        Ne[j] <- N[j]-NO[j]
      }
    }
    
    er2 <- array(NO,dim=c(48,7,8))
    Sij <- matrix(0,48,7)
    for (pp in 1:8){
      Sij <- Sij + er2[,,pp]
    }
    
    Sj <- colSums(Sij)
    S <- sum(Sj)
    print(S)
    # calculate new lambda0
    aLS <- aLS + S
    bLT <- bLT + T 
    lambda0 <- rgamma(1, aLS, bLT)
    # calculate new lambdat
    alphadsj <- alphadsj + Sj 
    # calculate new lambdat
    delta <- rdirichlet(1, alphadsj) *7
    eta <- matrix(0,7,48)
    for (j in 1:7){
      alphansij[j,] <- alphansij[j,] + Sij[,j]
      eta[j,] <- rdirichlet(1, alphansij[j,])*48
    }
    
    lambdat <- lambda0  * eta * as.vector(delta)
    er <-  array(rep(0,48*7*5), dim = c(48,7,5))
    er[,,c(1:5)] <- t(lambdat)
    lambdatv <- as.vector(er)
    lambda0v <- c(lambda0v, lambda0)
    zv <- cbind(zv, c(A0[1,1],A0[2,2],A0[3,3]))
    print(lambda0)
    
    zt2 <- c(5,zt[-1680])
    bigZ <- matrix(0,3,3)
    for (p in 1:3){
      for (q in 1:3){
        bigZ[p,q] <- length(which(zt == (p-1) & zt2 == (q-1)))
      }
    }
    
    abcdZ1 <- abcdZ1 + bigZ[,1]
    abcdZ2 <- abcdZ2 + bigZ[,2]
    abcdZ3 <- abcdZ3 + bigZ[,3]
    
    A0[1,] <- rdirichlet(1,abcdZ1)
    A0[2,] <- rdirichlet(1,abcdZ2)
    A0[3,] <- rdirichlet(1,abcdZ3)
    print(A0)
    print(kk)
  }
  
  List <- list("lv"= lambda0v, "zv"=zv, "z"=zt, "ltv" = lambdatv, "a0" = A0)
  return(List)
  
}
obs <- sort(N, decreasing = FALSE) # c(rep(0,64), rep(1,17), rep(2,10), rep(3,6), rep(4,3))
dzip <- function (x, mu, sigma) {
  ifelse((x == 0), (sigma + (1 - sigma) * exp(-mu)), ((1 - sigma) * dpois(x,mu)))
}
fit_zip = fitdistr(obs, dzip, start = list(mu = 60, sigma = 0.99), lower = list(p = 0.00001))
sigma <- as.numeric(fit_zip$estimate[2])
cc <- MMPP3(N, aL = sum(N), bL = 1680, sig = sigma)


plot(N[c(1:1344)], pch = ".", xlab = "",ylim = c(0,max(N)*2), ylab = "")
lines(N[c(1:1344)])
par(new=TRUE)
plot(b$z[c(1:1344)], pch = '.', col = "red", ylab = "", xlab = "",
     ylim = c(-1,2.2), yaxt = "n") #, main = "C18554")
lines(b$z[c(1:1344)], col = 'red')
par(new = TRUE)
plot(b$ltv[c(1:1344)], pch = '.', col = 'green', ylab = '', ylim = c(0,max(N)*2),
     xlab = '', yaxt = 'n')
lines(b$ltv[c(1:1344)], col = "green")


plot(N[c(1345:1680)], pch = ".", xlab = "",ylim = c(0,max(N)*2), ylab = "")
lines(N[c(1345:1680)])
par(new=TRUE)
plot(b$z[c(1345:1680)], pch = '.', col = "red", ylab = "", xlab = "",
     ylim = c(-1,2.2), yaxt = "n") #, main = "C18554")
lines(b$z[c(1345:1680)], col = 'red')
par(new = TRUE)
plot(b$ltv[c(1345:1680)], pch = '.', col = 'green', ylab = '', ylim = c(0,max(N)*2),
     xlab = '', yaxt = 'n')
lines(b$ltv[c(1345:1680)], col = "green")




plot(N, pch = ".", xlab = "Time bins",ylim = c(0,max(N)*2), ylab = "")
lines(N)
par(new=TRUE)
plot(c$z, pch = '.', col = "red", ylab = "", xlab = "",
     ylim = c(-3,3), yaxt = "n") #, main = "C18554")
lines(c$z, col = 'red')
par(new = TRUE)
plot(b$ltv, pch = '.', col = 'green', ylab = '', ylim = c(0,max(N)*2),
     xlab = '', yaxt = 'n')
lines(b$ltv, col = "green")



plot(N[c(1:336)], type = 'l', xlab = "",ylim = c(0,max(N)*2), ylab = "")
par(new=TRUE)
plot(c$z[c(1:336)], type = 'l', col = "red", ylab = "", xlab = "",
     ylim = c(-5,5), yaxt = "n") #, main = "C18554")
par(new = TRUE)
plot(c$ltv[c(1:336)], type = 'l', col = 'green', ylab = '', ylim = c(0,max(N)*2),
     xlab = '', yaxt = 'n')


