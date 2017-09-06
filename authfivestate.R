aZ0 <- 30; bZ0 <- 30; cZ0 <- 30; dZ0 <- 30; eZ0 <- 30
aZ1 <- 30; bZ1 <- 30; cZ1 <- 30; dZ1 <- 30; eZ1 <- 30
aZ2 <- 30; bZ2 <- 30; cZ2 <- 30; dZ2 <- 30; eZ2 <- 30
aZ3 <- 30; bZ3 <- 30; cZ3 <- 30; dZ3 <- 30; eZ3 <- 30
aZ4 <- 30; bZ4 <- 30; cZ4 <- 30; dZ4 <- 30; eZ4 <- 30

MMPP5 <- function(N, aL, bL,sig){
  q90 <- as.numeric(quantile(N, probs = 0.60))
  q99 <- as.numeric(quantile(N, probs = 0.99))
  aE1 <- q90
  aE2 <- q99
  bE1 <- 2
  bE2 <- 2
  # dirichlet for day of the week
  alphad <- matrix(sum(N)/7,1,7)
  # dirichlet for time of the day
  alphan <- matrix(sum(N)/336,7,48)
  # set initial value for lambda0
  lambda0 <- mean(N) #rgamma(1,aL, bL)
  set.seed(204)
  # day of week effect
  delta <- rdirichlet(1,alphad)*7
  set.seed(205)
  # time of day effect
  eta <- matrix(0,7,48)
  for (j in 1:7){
    eta[j,] <- rdirichlet(1, alphan[j,] )*48
  }
  # compute the periodic intensities
  # organises lambda_t into a 7*48 matrix 
  # row is day of the week and column is time of the day
  lambdat <- lambda0 * eta * as.vector(delta)
  # replicate this matrix eight times 
  lambdatoarray <-  array(rep(0,48*7*8), dim = c(48,7,8))
  lambdatoarray[,,c(1:8)] <- t(lambdat)
  lambdatv <- as.vector(lambdatoarray)
  # set initial value for transition matrix
  A0 <- rdirichlet(5,c(1/5,1/5,1/5,1/5,1/5))
  # compute the mean number of events for each day
  lambdahat <- as.numeric(0)
  for (dd in 1:56){
    lambdahat[dd] <- mean(N[c((48*(dd-1)+1):(48*dd))])
  }
  # computes which day each time bin belongs to
  days <- (c(1:2688)-1)%/%48 + 1
  aLS <- aL ; bLT <- bL
  abcdZ1 <- matrix(2688/25,1,5); abcdZ2 <- matrix(2688/25,1,5)
  abcdZ3 <- matrix(2688/25,1,5); abcdZ4 <- matrix(2688/25,1,5)
  alphadsj <- alphad
  alphansij <- alphan
  
  # store the lambda0 and transition probabilities for 
  # each iteration
  lambda0v <- c()
  zv <- cbind()
  # iteration begins
  for (kk in 1:150) {
    # function that calculates poisson-negative binomial
    poisnbin <- function(i,tt,aE,bE) {
      return(dpois((N[tt]-i), lambdatv[tt])*dnbinom(i, aE, bE/(1+bE)))}
    # function that sums up the probabilities of all
    # possible number of periodic events
    poisbinl <- function(tt,aE,bE){
      b <- as.numeric(0)
      for (j in 0:N[tt]){
        b[j+1] <- poisnbin(j,tt,aE,bE)
      }
      return(sum(b))
    }
    
    # likelihood function
    Pntzt <- function(tt,z){
      likelihood <- 0
      if (N[tt]==0 & lambdahat[days[tt]]==0){
        if (z == 0){likelihood <- 1}
        if (z == 1){likelihood <- 0}
        if (z == 2){likelihood <- 0}
        if (z == 3){likelihood <- 0}
        if (z == 4){likelihood <- 0}} 
      else if (N[tt]==0 & lambdahat[days[tt]]> 0){
        if (z == 0){likelihood <- sig}
        if (z == 1){likelihood <- (1-sig)*exp(-lambdahat[days[tt]])}
        if (z == 2){likelihood <- (1-sig)*exp(-lambdatv[tt])}
        if (z == 3){likelihood <- 0}
        if (z == 4){likelihood <- 0} }
      else {
        if (z == 0){likelihood <- 0}
        if (z == 1){likelihood <- dpois(N[tt], lambdahat[days[tt]])}
        if (z == 2){likelihood <- dpois(N[tt], lambdatv[tt])}
        if (z == 3){likelihood <- poisbinl(tt,aE1,bE1)}
        if (z == 4){likelihood <- poisbinl(tt,aE2,bE2)} }
      return(likelihood)
    }
    # Define the transtion probabilities 
    pznznm10 <- A0[,1]
    pznznm11 <- A0[,2]
    pznznm12 <- A0[,3]
    pznznm13 <- A0[,4]
    pznznm14 <- A0[,5]
    # Forward recursion
    alphazn1v<- matrix(0,5,2688)
    # Set initial value of forward recursion 
    initialalpha <- c(1/5,1/5,1/5,1/5,1/5) * 
      c(Pntzt(1,0),Pntzt(1,1),Pntzt(1,2),Pntzt(1,3),Pntzt(1,4))
    alphazn1v[,1] <- initialalpha/sum(initialalpha)
    for (t in 2:2688){
      if (N[t]==0 & lambdahat[days[t]] == 0){
        a1 <- sum(alphazn1v[,(t-1)]* pznznm10)
        a2 <- 0; a3 <- 0; a4 <- 0 ; a5 <- 0
      }
      else if (N[t]==0 & lambdahat[days[t]]> 0){
        a1 <- sig* sum(alphazn1v[,(t-1)]* pznznm10)
        a2 <- (1-sig)*exp(-lambdahat[days[t]])*sum(alphazn1v[,(t-1)]* pznznm11)
        a3 <- (1-sig)*exp(-lambdatv[t])*sum(alphazn1v[,(t-1)]* pznznm12)
        a4 <- 0 ; a5 <- 0
      } else {
        a1 <- 0
        a2 <- dpois(N[t], lambdahat[days[t]]) * sum(alphazn1v[,(t-1)]* pznznm11)
        a3 <- dpois(N[t], lambdatv[t]) * sum(alphazn1v[,(t-1)]* pznznm12)
        a4 <- poisbinl(t,aE1,bE1) * sum(alphazn1v[,(t-1)]* pznznm13)
        a5 <- poisbinl(t,aE2,bE2) * sum(alphazn1v[,(t-1)]* pznznm14)
      }
      alphazn1v[1,t] <- a1/(a1+a2+a3+a4+a5)
      alphazn1v[2,t] <- a2/(a1+a2+a3+a4+a5)
      alphazn1v[3,t] <- a3/(a1+a2+a3+a4+a5)
      alphazn1v[4,t] <- a4/(a1+a2+a3+a4+a5) 
      alphazn1v[5,t] <- a5/(a1+a2+a3+a4+a5) 
    }
    
    # backward recursion
    # set initial value of backward recursion
    betazN <- 1/5
    betazn1v <- matrix(0,5,2688)
    betazn1v[,2688] <- betazN
    for (t in 2687:1){
      if (N[t+1]==0 & lambdahat[days[t+1]] == 0) {
        b <- as.numeric(0)
        for (g in 1:5){
          b[g] <- betazn1v[1,t+1] * A0[g,1]
        }
      } else if (N[t+1]==0 & lambdahat[days[t]] > 0) {
        b <- as.numeric(0)
        for (g in 1:5){
          b[g] <- betazn1v[1,t+1] * A0[g,1] +
            betazn1v[2,t+1] * exp(-lambdahat[days[t+1]]) * A0[g,2] + 
            betazn1v[3,t+1] * exp(-lambdatv[t+1]) * A0[g,3]
        }
      } else {
        b <- as.numeric(0)
        for (g in 1:5) {
          b[g] <- betazn1v[2,t+1] * dpois(N[t+1], lambdahat[days[t+1]]) * A0[g,2]+
          betazn1v[3,t+1] * dpois(N[t+1], lambdatv[t+1]) * A0[g,3] +
            betazn1v[4,t+1] * poisbinl(t+1,aE1,bE1) * A0[g,4] +
            betazn1v[5,t+1] * poisbinl(t+1,aE2,bE2) * A0[g,5]
        }
      }
      betazn1v[1,t] <- b[1]/sum(b)
      betazn1v[2,t] <- b[2]/sum(b)
      betazn1v[3,t] <- b[3]/sum(b)
      betazn1v[4,t] <- b[4]/sum(b)
      betazn1v[5,t] <- b[5]/sum(b)
    }
    
    gamma <- matrix(0,5,2688)
    g1 <- alphazn1v[1,]*betazn1v[1,] 
    g2 <- alphazn1v[2,]*betazn1v[2,] 
    g3 <- alphazn1v[3,]*betazn1v[3,]
    g4 <- alphazn1v[4,]*betazn1v[4,]
    g5 <- alphazn1v[5,]*betazn1v[5,]
    gamma[1,] <- g1 / (g1+g2+g3+g4+g5)
    gamma[2,] <- g2 / (g1+g2+g3+g4+g5)
    gamma[3,] <- g3 / (g1+g2+g3+g4+g5)
    gamma[4,] <- g4 / (g1+g2+g3+g4+g5)
    gamma[5,] <- g5 / (g1+g2+g3+g4+g5)
    
    pztzt1 <- array(0,c(5,5,2688))
    pztzt1v <- matrix(0,5,5)
    for (t in 2687:1){
      for (p in 1:5){
        for (q in 1:5){
          pztzt1v[p,q] <- alphazn1v[p,t]  * Pntzt(t+1,(q-1)) *
            A0[q,p] *betazn1v[q,t+1]
        }
      }
      
      p1 <- sum(pztzt1v[,1]) # p(zt1 = 0)
      p2 <- sum(pztzt1v[,2]) # p(zt1 = 1)
      p3 <- sum(pztzt1v[,3]) # p(zt1 = 2)
      p4 <- sum(pztzt1v[,4])
      p5 <- sum(pztzt1v[,5])
      if (p1 == 0) {pztzt1v[,1] <- A0[1,];p1 = 1}
      if (p2 == 0) {pztzt1v[,2] <- A0[2,];p2 = 1}
      if (p3 == 0) {pztzt1v[,3] <- A0[3,];p3 = 1}
      if (p4 == 0) {pztzt1v[,4] <- A0[4,];p4 = 1}
      if (p5 == 0) {pztzt1v[,5] <- A0[5,];p5 = 1}
      p <- c(p1,p2,p3,p4,p5)
      for (q in 1:5){
        pztzt1[,q,t] <- pztzt1v[,q]/p[q]
        if (sum(pztzt1[,q,t])==0){pztzt1[,q,t] = 1/5}
      }
    }
    #print(pztzt1[,,c(1:10)])
    zt <- as.numeric(0)
    zt[2688] <- sample(c(0,1,2,3,4),1, prob = gamma[,2688])
    for (tt in 2687:1){
      if (zt[tt+1] ==0 ){zt[tt] <- sample(c(0,1,2,3,4), 1,prob = pztzt1[,1,tt])}
      else if (zt[tt+1] ==1 ){zt[tt] <- sample(c(0,1,2,3,4), 1,prob = pztzt1[,2,tt])}
      else if (zt[tt+1] ==2) {zt[tt] <- sample(c(0,1,2,3,4), 1,prob = pztzt1[,3,tt])}
      else if (zt[tt+1] ==3) {zt[tt] <- sample(c(0,1,2,3,4), 1,prob = pztzt1[,4,tt])}
      else {zt[tt] <- sample(c(0,1,2,3,4), 1, prob = pztzt1[,5,tt])}
    }
    print(table(zt))
    
    
    samplePNB <- function(tt,aE,bE){
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
    Np <- as.numeric(0)
    for (j in 1:2688){
      if (zt[j]==0){
        NO[j] <- 0
        Ne[j] <- 0
        Np[j] <- 0
      }
      else if (zt[j]==1){
        Np[j] <- N[j]
        Ne[j] <- 0
        NO[j] <- 0
      } 
      else if (zt[j]==2){
        NO[j] <- N[j]
        Ne[j] <- 0
        Np[j] <- 0
      }
      else if (zt[j] == 3){
        NO[j] <- samplePNB(j,aE1,bE1)
        Ne[j] <- N[j]-NO[j]
        Np[j] <- 0
      }
      else if (zt[j] == 4){
        NO[j] <- samplePNB(j,aE2,bE2)
        Ne[j] <- N[j]-NO[j]
        Np[j] <- 0
      }
      
    }
    
    NOtoarray <- array(NO,dim=c(48,7,8))
    # count the number of periodic events for each time bin
    Sij <- matrix(0,48,7)
    for (pp in 1:8){
      Sij <- Sij + NOtoarray[,,pp]
    }
    # marginalise out all values for weeks and days
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
    lambdatoarray <-  array(rep(0,48*7*8), dim = c(48,7,8))
    lambdatoarray[,,c(1:8)] <- t(lambdat)
    lambdatv <- as.vector(lambdatoarray)
    #lambdatv[which(zt==1)] <- lambdahat[days[which(zt==1)]]
    lambda0v <- c(lambda0v, lambda0)
    zv <- cbind(zv, c(A0[1,1],A0[2,2],A0[3,3],A0[4,4]))
    print(lambda0)
    # calculate the number of the cases zt|zt+1
    zt2 <- c(5,zt[-2688])
    bigZ <- matrix(0,5,5)
    for (p in 1:5){
      for (q in 1:5){
        bigZ[p,q] <- length(which(zt == (p-1) & zt2 == (q-1)))
      }
    }
    # generate new transition probabilities and replace the matrix
    
    abcdZ1 <- abcdZ1 + bigZ[,1]
    abcdZ2 <- abcdZ2 + bigZ[,2]
    abcdZ3 <- abcdZ3 + bigZ[,3]
    abcdZ4 <- abcdZ4 + bigZ[,4]
    abcdZ5 <- abcdZ4 + bigZ[,5]
    
    A0[1,] <- rdirichlet(1,abcdZ1)
    A0[2,] <- rdirichlet(1,abcdZ2)
    A0[3,] <- rdirichlet(1,abcdZ3)
    A0[4,] <- rdirichlet(1,abcdZ4)
    A0[5,] <- rdirichlet(1,abcdZ5)
    
    print(A0)
  }
  
  List <- list("lv"= lambda0v, "zv"=zv, "z"=zt, "ltv"= lambdatv, "a0" = A0)
  return(List)
}



obs <- sort(N, decreasing = FALSE) # c(rep(0,64), rep(1,17), rep(2,10), rep(3,6), rep(4,3))
dzip <- function (x, mu, sigma) {
  ifelse((x == 0), (sigma + (1 - sigma) * exp(-mu)), ((1 - sigma) * dpois(x,mu)))
}

fit_zip = fitdistr(obs, dzip, start = list(mu = 30, sigma = 0.9), lower = list(p = 0.00001))
sigma <- as.numeric(fit_zip$estimate[2])
ee <- MMPP5(N, aL = sum(N), bL = 2688, sig = sigma)


   par(mfrow = c(1,1))
plot(N, pch = ".", xlab = "", ylim =c(0,max(N))*2)
lines(N)
par(new=TRUE)
plot(ee$z, pch = '.', col = "red", ylab = "", xlab = "",
     ylim = c(-10,5), yaxt = "n", main = "C10123 C801, aE = 50, bE = 2")
lines(ee$z, col = 'red')



