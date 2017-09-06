library("Matrix", lib.loc="C:/Program Files/R/R-3.3.1/library")
library(MCMCpack)

### fixed starting values
# gamma distribution for lambda0
aL <- 50000
bL <- 3000
# dirichlet for day of the week
alphad <- matrix(sum(N)/7,1,7)
# dirichlet for time of the day
alphan <- matrix(sum(N)/336,7,48)
# beta distribution for z0,z1
aZ0 <- 100
bZ0 <- 20000
aZ1 <- 20000
bZ1 <- 100
# gamma distrubution for poisson distribution for events
aE <- 2
bE <- 1
# total number of time slots and timeslots per day
T <- 2688
D=48



########### Now run algorithm ############

MMPP <- function(N, aL, bL, sig){
  aE <- as.numeric(quantile(N, probs = 0.9))
  bE <- 2
  ## initial values ##
  # initial distribution
  pi0 <- c(0.5,0.5)
  # initial lambdat 
  lambda0 <- mean(N) #rgamma(1,aL, bL)
  # dirichlet for day of the week
  alphad <- matrix(sum(N)/7,1,7)
  # dirichlet for time of the day
  alphan <- matrix(sum(N)/336,7,48)
  set.seed(204)
  delta <- rdirichlet(1,alphad)*7
  set.seed(205)
  eta <- matrix(0,7,48)
  for (j in 1:7){
    eta[j,] <- rdirichlet(1, alphan[j,] )*48
  }
  lambdat <- lambda0 * eta * as.vector(delta)
  er <-  array(rep(0,48*7*8), dim = c(48,7,8))
  er[,,c(1:8)] <- t(lambdat)
  lambdatv <- as.vector(er)
  
  # initial transition matrix for hidden markov chain
  set.seed(302)
  A0 <- matrix(1/2,2,2)
  aLS <- aL ; bLT <- bL
  abcdZ1 <- matrix(2688/4,1,2); abcdZ2 <- matrix(2688/4,1,2)
  abcdZ3 <- matrix(2688/4,1,2); abcdZ4 <- matrix(2688/4,1,2)
  alphadsj <- alphad
  alphansij <- alphan
  
  ## loop ##
  
  lambda0v <- c()
  zv <- cbind()
  for (kk in 1:150){
    ## likelihood function
    PNtzt <- function(t,z){
      likelihood <- as.numeric(0)
      if (z == 0 ){
        if (N[t] ==0) {
          likelihood <- sig + (1-sig)*exp(-lambdatv[t])
        } else {
          likelihood <-(1-sig)*dpois(N[t],lambdatv[t])
        }
      } else if (z == 1 ){
        b <- as.numeric(0)
        c <- as.numeric(0)
        for (i in 0:N[t]){
          b[i+1] <- dpois((N[t]-i),lambdatv[t])
          c[i+1] <- dnbinom(i, aE,bE/(1+bE))
        }
        likelihood <- sum(b*c)
      }
      return(likelihood)
    } 
    # first value for forward recursion
    alphaz1<- pi0*c(PNtzt(1,0), PNtzt(1,1))/(pi0[1]*PNtzt(1,0) + pi0[2]* PNtzt(1,1))
    # first value for backward recursion
    betazN <- 0.5
    # probabilities for landing in each state given previous state
    pznznm10 <- A0[,1]
    pznznm11 <- A0[,2]
    # begins forward backward algorithm
    alphazn1v <-matrix(0,2,2688)
    alphazn1v[,1] <- alphaz1
    for (t in 2:2688){
      a1 <- PNtzt(t,0) *sum(alphazn1v[,(t-1)] *pznznm10)
      a2 <- PNtzt(t,1) * sum(alphazn1v[,(t-1)] * pznznm11)
      alphazn1v[1,t] <- a1/(a1+a2)
      alphazn1v[2,t] <- a2/(a1+a2)
    }
    #print(alphazn1v[,20])
    #print(alphazn1v[1,c(1:5)])
    betazn1v <- matrix(0,2,2688)
    betazn1v[,2688] <- betazN
    for (t in 2687:1){
      b1 <- betazn1v[1,t+1]*PNtzt(t+1,0)*A0[1,1] +
        betazn1v[2,t+1]*PNtzt(t+1,1)* A0[1,2]
      b2 <- betazn1v[1,t+1]*PNtzt(t+1,0)*A0[2,1] +
        betazn1v[2,t+1]*PNtzt(t+1,1)*A0[2,2]
      betazn1v[1,t] <- b1 / (b1+b2)
      betazn1v[2,t] <- b2/ (b1+b2)
    }
    #pN <- colSums(alphazn1v*betazn1v)
    gamma <- matrix(0,2,2688)
    g1 <- alphazn1v[1,]*betazn1v[1,] #/ pN
    g2 <- alphazn1v[2,]*betazn1v[2,] #/ pN
    gamma[1,] <- g1 / (g1+g2)
    gamma[2,] <- g2 / (g1+g2)
    
    ####### p(zt|zt+1, N) ######
    pztzt1 <- matrix(0,4,2687)
    for (t in 1:2687){
      pzt0zt11 <- alphazn1v[1,t]  * PNtzt(t+1,1) * A0[1,2] *betazn1v[2,t+1]
      pzt1zt11 <- alphazn1v[2,t]  * PNtzt(t+1,1) * A0[2,2] *betazn1v[2,t+1]
      pzt0zt10 <- alphazn1v[1,t]  * PNtzt(t+1,0) * A0[1,1] *betazn1v[1,t+1]
      pzt1zt10 <- alphazn1v[2,t]  * PNtzt(t+1,0) * A0[2,1] *betazn1v[1,t+1]
      pztzt1all <- pzt0zt11 + pzt1zt11 + pzt0zt10 + pzt1zt10
      p1 <- pzt0zt10/pztzt1all # p(zt = 0, zt1 = 0)
      p2 <- pzt1zt10/pztzt1all # p(zt = 1, zt1 = 0)
      p3 <- pzt0zt11/pztzt1all # p(zt = 0, zt1 = 1)
      p4 <- pzt1zt11/pztzt1all # p(zt = 1, zt1 = 1)
      pztzt1[1,t] <- p1/(p1+p2) # p(zt = 0 | zt1 = 0)
      pztzt1[2,t] <- p2/(p1+p2) # p(zt = 1 | zt1 = 0)
      pztzt1[3,t] <- p3/(p3+p4) # p(zt = 0 | zt1 = 1)
      pztzt1[4,t] <- p4/(p3+p4) # p(zt = 1 | zt1 = 1)
      
     # if (alphazn1v[1,t] >0 && alphazn1v[2,t] >0 ){
    }
    #print(pztzt1[,c(1:5)])
    # generate from p(zt|zt+1, N)
    zt <- as.numeric(0)
    zt[2688] <- sample(c(0,1),1, prob= gamma[,2688])
    for (tt in 2687:1){
      if (zt[tt+1] ==0 ){zt[tt] <- sample(c(0,1), 1,prob = pztzt1[c(1,2),tt])}
      else if (zt[tt+1] ==1 ){zt[tt] <- sample(c(0,1), 1,prob = pztzt1[c(3,4),tt])}
    }
    print(table(zt))
    # function that samples number of periodic events
    samplePNB <- function(tt){
      if (N[tt]>0){
        poisnbin <- function(i){
          return(dpois(N[tt] - i, lambdatv[tt])* dnbinom(i,aE, bE/(1+bE)))}
        pbv <- as.numeric(0)
        for (i in 0:(N[tt])){
          pbv[i+1] <- poisnbin(i)}
        probpbv <- pbv/sum(pbv)
        #plot(probpbv, pch =19)
        return(sample(c(0:N[tt]), 1, prob = probpbv))}
      else {return(0)}}
          # sample the number of periodic events and rare events 
        NO <- as.numeric(0)
        Ne <- as.numeric(0)
        for (j in 1:2688){
          if (zt[j]==0){
            NO[j] <- N[j]
            Ne[j] <- 0}
          else if (zt[j] == 1){
            NO[j] <- samplePNB(j)
            Ne[j] <- N[j]-NO[j]}
          }
    # calculates the number of distributions
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
    
    # now calculate the transition probabilities
    zt2 <- c(5,zt[-2688])
    bigZ <- matrix(0,2,2)
    for (p in 1:2){
      for (q in 1:2){
        bigZ[p,q] <- length(which(zt == (p-1) & zt2 == (q-1)))
      }
    }
    abcdZ1 <- abcdZ1 + bigZ[,1]
    abcdZ2 <- abcdZ2 + bigZ[,2]
    
    z0 <- rbeta(1,abcdZ1[2],abcdZ1[1])
    z1 <- rbeta(1,abcdZ2[1],abcdZ2[2])
    
    ## replace old values with new values
    
    A0 <- matrix(c(1-z0, z1, z0, 1-z1),2,2)
    
    er <-  array(rep(0,48*7*8), dim = c(48,7,8))
    er[,,c(1:8)] <- t(lambdat)
    lambdatv <- as.vector(er)
    lambda0v <- c(lambda0v, lambda0)
    zv <- cbind(zv, c(z0,z1))
    print(A0)
    print(lambda0)
  }
  List <- list("lv"= lambda0v, "zv"=zv, "z"=zt, "ltv" = lambdatv, 'a0' = A0 )
  return(List)
}


obs <- sort(N, decreasing = FALSE) # c(rep(0,64), rep(1,17), rep(2,10), rep(3,6), rep(4,3))
dzip <- function (x, mu, sigma) {
  ifelse((x == 0), (sigma + (1 - sigma) * exp(-mu)), ((1 - sigma) * dpois(x,mu)))
}
fit_zip = fitdistr(obs, dzip, start = list(mu = 30, sigma = 0.9), lower = list(p = 0.00001))
sigma <- as.numeric(fit_zip$estimate[2])


bb <- MMPP(N, aL= sum(N), bL= 2688, sig = sigma )

par(mfrow = c(2,2))
for (jj in 1:8){
  plot(a$ltv[c(1:336)], pch = '.', ylim = c(-0.01,max(N)), ylab = "")
  lines(a$ltv[c(1:336)])
  par(new = TRUE )
  plot(N[c((336*(jj-1)+1):(jj*336))], pch = '.', col = 'red', ylim = c(-0.01,max(N)),
       main = paste("week",jj), ylab = "")
  lines(N[c((336*(jj-1)+1):(jj*336))], col = "red")
  par(new = TRUE)
  plot(a$NE[c((336*(jj-1)+1):(jj*336))], pch = '.', ylim = c(-0.01,max(N)), 
       yaxt = 'n', ylab = "", col = 'blue')
  lines(a$NE[c((336*(jj-1)+1):(jj*336))], col = 'blue')
  par(new = TRUE)
  plot(a$NO[c((336*(jj-1)+1):(jj*336))], pch = '.', ylim = c(-0.01,max(N)), 
       yaxt = 'n', ylab = "", col = 'purple')
  lines(a$NO[c((336*(jj-1)+1):(jj*336))], col = 'purple')
  par(new = TRUE)
  plot(a$z[c((336*(jj-1)+1):(jj*336))], col = "green", pch = '.', yaxt = 'n',
       ylim = c(-2,1.5), ylab ='')
  lines(a$z[c((336*(jj-1)+1):(jj*336))], col = "green")
}






par(mfrow = c(1,1))
plot(a$lv, pch = ".", ylab = "Simulated lambda0", xlab = "Iteration", 
     main = "lambda0 = 2, z0 = 0.4, z1 = 0.9", ylim = c(min(a$lv), max(a$lv)))
lines(a$lv)
abline(h = 10, col = "red")
plot(a$zv[1,], pch = ".", ylim = c(0,1), col = "green", ylab = "",
     main = "lambda0 = 2, z0 = 0.4, z1 = 0.9")
lines(a$zv[1,], col = "green")
abline(h=0.3, col = "green" )
par(new = TRUE)
plot(a$zv[2,], pch = ".", ylim = c(0,1), col = "red", ylab = "")
lines(a$zv[2,], col = "red")
abline(h=0.8, col = "red")
par(mfrow = c(2,1))
par(mfrow = c(1,1))
plot(N[c(1:1344)], pch = ".", ylim = c(0,max(N)*2), xlab = "Time bins", ylab = '') #  main = "C1085 and C553",
lines(N[c(1:1344)])
par(new = TRUE)
plot(a$z[c(1:1344)], pch = ".", col = "red", ylab = "", xlab = "",
     ylim = c(-2,2), yaxt = "n")
lines(a$z[c(1:1344)], col = "red")
par(new = TRUE)
plot(a$ltv[c(1:1344)], pch = '.', col = 'green', ylab = '', ylim = c(0,max(N)*2),
     xlab = '')
lines(a$ltv[c(1:1344)], col = "green")


par(mfrow = c(1,1))
plot(N[c(1345:2688)], pch = ".", ylim = c(0,max(N)*2), xlab = "Time bins", ylab = '') #  main = "C1085 and C553",
lines(N[c(1345:2688)])
par(new = TRUE)
plot(a$z[c(1345:2688)], pch = ".", col = "red", ylab = "", xlab = "",
     ylim = c(-2,2), yaxt = "n")
lines(a$z[c(1345:2688)], col = "red")
par(new = TRUE)
plot(a$ltv[c(1345:2688)], pch = '.', col = 'green', ylab = '', ylim = c(0,max(N)*2),
     xlab = '')
lines(a$ltv[c(1345:2688)], col = "green")







par(mfrow = c(2,1))
plot(z[c(1:100)], pch = ".")
lines(z)
plot(a$z[c(1:100)], col = "red", pch = '.')


plot(a$ltv[c(1:336)], pch = '.', ylim = c(0,50))
lines(a$ltv[c(1:336)])
par(new = TRUE )
plot(N[c(673:1008)], pch = '.', col = 'red', ylim = c(0,50))
lines(N[c(673:1008)], col = "red")

par(new = FALSE)
par(mfrow = c(1,1))
plot(a$NO[c(1:336)], pch = '.', ylim = c(-0.01,max(N)), col = 'red')
lines(a$NO[c(1:336)], col = 'red')
par(new = TRUE)
plot(N[c(1:336)], pch = '.', ylim = c(-0.01,max(N)), col = 'green')
lines(N[c(1:336)], col = 'green')
par(new = TRUE)
plot(a$NE[c(1:336)], pch = '.', ylim = c(-0.01,max(N)), col = 'blue')
lines(a$NE[c(1:336)], col = 'blue')

