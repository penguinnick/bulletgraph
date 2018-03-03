#make matrix
m <- matrix(c(162,49,57,40,43,49),2,3,byrow = T)
#get proportions
m.prop <- prop.table(m,1)
datG=m[1,] #zir.gray$Zircon
datR=m[2,] #zir.black$Zircon
dat <- m.prop

bulletplot=function(datG,datR){
  #m <- matrix(c(162,49,57,40,43,49),2,3,byrow = T)
  perrorRange <- function(a, b, c){   #-- Return Standard Error. Let a=count or proportion, b=total, c= confidence level.
    #error Range function for proportions of population: 
    p <- a/b   #proportion
    q <- 1-p   #remaining proportion
    se <- sqrt(p*q)/sqrt(b) # standard error
    t <-  qt(c, (b-1))
    t*se
  }
  #get quantiles
  #qz=quantile(dat,p=c(0.005,0.025,0.10,0.90,0.975,0.995))
  #calculate confidence intervals
  #t-values for 80, 95, 99% intervals
  g <- sum(datG)
  r <- sum(datR)
  #proportions at type for site g and proportions at type for site r
  g.p <- datG/g
  r.p <- datR/r
  p.1 <- c(g.p,r.p)
  #error ranges 80%, 95%, 99%
  g.mod <- sapply(1:3, function(x,y) perrorRange(y[1,x],g,.90),y=m)
  g.hi  <- sapply(1:3, function(x,y) perrorRange(y[1,x],g,.975),y=m)
  g.exhi<- sapply(1:3, function(x,y) perrorRange(y[1,x],g,.995),y=m)
  g.er <- matrix(c(g.mod,g.hi,g.exhi),3,3,byrow = T)
  g.max80 <- g.p + g.er[1,] 
  g.min80 <- g.p - g.er[1,]
  g.max95 <- g.p + g.er[2,]
  g.min95 <- g.p - g.er[2,]
  g.max99 <- g.p + g.er[3,]
  g.min99 <- g.p - g.er[3,]
  
  r.mod <- sapply(1:3, function(x,y) perrorRange(y[2,x],r,.90),y=m)
  r.hi  <- sapply(1:3, function(x,y) perrorRange(y[2,x],r,.975),y=m)
  r.exhi<- sapply(1:3, function(x,y) perrorRange(y[2,x],r,.995),y=m)
  r.er <- matrix(c(r.mod,r.hi,r.exhi),3,3,byrow = T)
  r.max80 <- r.p + r.er[1,] 
  r.min80 <- r.p - r.er[1,]
  r.max95 <- r.p + r.er[2,]
  r.min95 <- r.p - r.er[2,]
  r.max99 <- r.p + r.er[3,]
  r.min99 <- r.p - r.er[3,]
  p.80min <- c(g.min80,r.min80)
  p.80max <- c(g.max80,r.max80)
  p.95min <- c(g.min95,r.min95)
  p.95max <- c(g.max95,r.max95)
  p.99min <- c(g.min99,r.min99)
  p.99max <- c(g.max99,r.max99)
  # #get range
  rz <- range(p.99min,p.99max)
  plot(1,mean(p.1),ylim=rz,pch=".",col='black',xlim=c(0,7),xlab="Gray is for Grander, Black is for Rawlins",ylab="Error Ranges")
  for (i in 1:7) {
    if ( i %% 1) {
      next
    }
    sl1 <- i-.2500
    sl2 <- i+.2500
    a <- i-.1250   
    b <- i+.1250
    c <- i-.0675
    d <- i+.0675
    e <- i-0.0225
    f <- i+0.0225
    if(i <= 3){
      j <- "gray"
    } else {
      j<-"black"
    }
    segments(sl1,p.1[i],sl2,p.1[i],col = j,lwd=2)
    rect(a, p.80min[i], b, p.80max[i], density = -1, col = j, border = NA)
    rect(c, p.95min[i], d, p.95max[i], density = -1, col = j, border = NA)
    rect(e, p.99min[i], f, p.99max[i], density = -1, col = j, border = NA)
  }
}
