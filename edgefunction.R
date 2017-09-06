setwd('Desktop/Thesis')
lanl <- read.table("try4file.txt", sep = ",")

####### different source and destination #######
########### order by number of events ##########

diffc <- lanl[which(as.matrix(lanl[,4])!=as.matrix(lanl[,5])),]
diffcfreq <- as.data.frame(table(diffc[,c(4:5)]))
diffcfreqorder <- diffcfreq[order(diffcfreq$Freq, decreasing = TRUE),]

##### function that returns the number of events in each time bin #####

whichdedge <- function(tb, ee){
  timebin <- tb
  ntimebin <- 86400/timebin
  
  singleedge <- lanl[which(is.element(as.vector(lanl[,4]), diffcfreqorder$V4[ee])==TRUE &
                             is.element(as.vector(lanl[,5]), diffcfreqorder$V5[ee])==TRUE),]
  
  edgeday <- singleedge[,1]%/%86400
  edgedaysec <- singleedge[,1]%%86400
  edgedaybelong <- cbind(edgeday+1, edgedaysec+1)
  
  edgedf <- data.frame(edgedaybelong)
  edge1 <- matrix(0,58,ntimebin)
  for (d in 1:58){
    if (is.element(d,edgedf$X1)==TRUE ){
      for (j in 1:ntimebin){
        a <- edgedf$X2[which(edgedf$X1==d)]
        edge1[d,j] <- length(which(a >((j-1)*timebin) & a <= (j*timebin) ))
      }
    }
  }
  
  testedgeweek <- matrix(0, 8, ntimebin*7 )
  par(mfrow = c(2,4))
  for (w in 1:8){
    testedgeweek[w,] <- as.vector(t(edge1[c(((w-1)*7+1):(7*w)),]))
    plot(testedgeweek[w,], ylab = "events", ylim = c(0, max(edge1)),
         main = paste("week",w) , pch = '.')
    lines(testedgeweek[w,])
  }
  par(new = FALSE)
  N <- testedgeweek
  return(N)
}

d <- whichdedge(600, 2)
N <- as.vector(t(d))
(propzero <- length(which(N == 0))/length(N))


####### similarly same source and destination #######
############# order by number of events #############

samec <- lanl[which(as.matrix(lanl[,4])==as.matrix(lanl[,5])),]
samecfreq <- as.data.frame(table(samec[,c(4:5)]))
samecfreqorder <- samecfreq[order(samecfreq$Freq, decreasing = TRUE),]

##### function that returns the number of events in each time bin #####

whichsedge <- function(tb, ee){
  timebin <- tb
  ntimebin <- 86400/timebin
  
  singleedge <- lanl[which(is.element(as.vector(lanl[,4]), samecfreqorder$V4[ee])==TRUE &
                             is.element(as.vector(lanl[,5]), samecfreqorder$V5[ee])==TRUE),]
  
  edgeday <- singleedge[,1]%/%86400
  edgedaysec <- singleedge[,1]%%86400
  edgedaybelong <- cbind(edgeday+1, edgedaysec+1)
  
  edgedf <- data.frame(edgedaybelong)
  edge1 <- matrix(0,58,ntimebin)
  for (d in 1:58){
    if (is.element(d,edgedf$X1)==TRUE ){
      for (j in 1:ntimebin){
        a <- edgedf$X2[which(edgedf$X1==d)]
        edge1[d,j] <- length(which(a >((j-1)*timebin) & a <= (j*timebin) ))
      }
    }
  }
  
  testedgeweek <- matrix(0, 8, ntimebin*7 )
  par(mfrow = c(2,4))
  for (w in 1:8){
    testedgeweek[w,] <- as.vector(t(edge1[c(((w-1)*7+1):(7*w)),]))
    #plot(testedgeweek[w,], ylab = "events", ylim = c(0, max(edge1)),
    #     main = paste("week",w) , pch = '.')
    #lines(testedgeweek[w,])
  }
  par(new = FALSE)
  N <- testedgeweek
  return(N)
}


a <- whichsedge(1800,4)
N <- as.vector(t(a))
(propzero <- length(which(N == 0))/length(N))



########## for netflow

netflow <- read.table("netflow2.txt", sep = ",")

netflowsourcedest <- as.data.frame(table(netflow[,c(2,3)]))
netflowsourcedesto <- netflowsourcedest[order(netflowsourcedest$Freq, decreasing = TRUE),]



whichnedge <- function(tb, ee){
  timebin <- tb
  ntimebin <- 86400/timebin
  
  singleedge <- netflow[which(is.element(as.vector(netflow[,2]), netflowsourcedesto$V2[ee])==TRUE &
                                is.element(as.vector(netflow[,3]), netflowsourcedesto$V3[ee])==TRUE),]
  
  edgeday <- singleedge[,1]%/%86400
  edgedaysec <- singleedge[,1]%%86400
  edgedaybelong <- cbind(edgeday+1, edgedaysec+1)
  
  edgedf <- data.frame(edgedaybelong)
  edge1 <- matrix(0,58,ntimebin)
  for (d in 1:58){
    if (is.element(d,edgedf$X1)==TRUE ){
      for (j in 1:ntimebin){
        a <- edgedf$X2[which(edgedf$X1==d)]
        edge1[d,j] <- length(which(a >((j-1)*timebin) & a <= (j*timebin) ))
      }
    }
  }
  
  testedgeweek <- matrix(0, 5, ntimebin*7 )
  par(mfrow = c(2,4))
  for (w in 1:5){
    testedgeweek[w,] <- as.vector(t(edge1[c(((w-1)*7+1):(7*w)),]))
    #plot(testedgeweek[w,], ylab = "events", ylim = c(0, max(edge1)),
    #   main = paste("week",w) , pch = '.')
    #lines(testedgeweek[w,])
  }
  par(new = FALSE)
  N <- testedgeweek
  return(N)
}

a <- whichnedge(1800,1)
N <- as.vector(t(a))
(propzero <- length(which(N == 0))/length(N))

par(mfrow = c(1,1))
plot(N, type = 'l')







