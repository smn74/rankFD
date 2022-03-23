########################################################################
#Function Effects
########################################################################


Effects<-function(response, factors,  effect = c("unweighted", 
                                                 "weighted")) 
{
  effect <- match.arg(effect)
  factorx <- factors
  samples <- split(response, factorx)
  fl <- levels(factorx)
  a <- nlevels(factorx)
  n <- sapply(samples,length)
  N <- sum(n)
  tmp1 <- sort(rep(1:a, a))
  tmp2 <- rep(1:a, a)
  
  p <- sapply(1:(a^2), function(arg) {
    x1 <- samples[[tmp1[arg]]]
    x2 <- samples[[tmp2[arg]]]
    rx1x2 <- rank(c(x1, x2))
    l1 <- length(x1)
    l2 <- length(x2)
    1/(l1 + l2) * (mean(rx1x2[(l1 + 1):(l1 + l2)]) - mean(rx1x2[1:l1])) + 
      0.5
  })
  
  pairRanks <- lapply(1:(a^2), function(arg) rank(c(samples[[tmp1[arg]]], 
        samples[[tmp2[arg]]])))
        

  
intRanks <- lapply(samples, rank)
  
placements <- lapply(1:(a^2), function(arg) 1/n[tmp1[arg]] * 
        (pairRanks[[arg]][(n[tmp1[arg]] + 1):(n[tmp1[arg]] + 
            n[tmp2[arg]])] - intRanks[[tmp2[arg]]]))

 
  
  V <- rep(0, a^4)
  help <- expand.grid(1:a, 1:a, 1:a, 1:a)
  h1 <- help[, 4]
  h2 <- help[, 3]
  h3 <- help[, 2]
  h4 <- help[, 1]
  for (u in 1:(a^4)) {
    i <- h1[u]
    j <- h2[u]
    r <- h3[u]
    s <- h4[u]
    if (i == r && j == s && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      ni <- length(xi)
      nj <- length(xj)
      ri <- rank(xi)
      rj <- rank(xj)
      rij <- rank(c(xi, xj))
      pj <- 1/ni * (rij[(ni + 1):(ni + nj)] - rj)
      pi <- 1/nj * (rij[1:ni] - ri)
      vi <- var(pi)/ni
      vj <- var(pj)/nj
      V[u] <- N * (vi + vj)
    }
    if (i == s && j == r && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      ni <- length(xi)
      nj <- length(xj)
      ri <- rank(xi)
      rj <- rank(xj)
      rij <- rank(c(xi, xj))
      pj <- 1/ni * (rij[(ni + 1):(ni + nj)] - rj)
      pi <- 1/nj * (rij[1:ni] - ri)
      vi <- var(pi)/ni
      vj <- var(pj)/nj
      V[u] <- -N * (vi + vj)
    }
    if (i == r && j != s && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      xs <- samples[[s]]
      ni <- length(xi)
      nj <- length(xj)
      ns <- length(xs)
      ri <- rank(xi)
      rj <- rank(xj)
      rs <- rank(xs)
      rij <- rank(c(xi, xj))
      ris <- rank(c(xi, xs))
      pij <- 1/nj * (rij[1:ni] - ri)
      pis <- 1/ns * (ris[1:ni] - ri)
      V[u] <- N * (cov(pij, pis)/ni)
    }
    if (i != r && j == s && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      xr <- samples[[r]]
      ni <- length(xi)
      nj <- length(xj)
      nr <- length(xr)
      ri <- rank(xi)
      rj <- rank(xj)
      rr <- rank(xr)
      rji <- rank(c(xj, xi))
      rjr <- rank(c(xj, xr))
      pji <- 1/ni * (rji[1:nj] - rj)
      prj <- 1/nr * (rjr[1:nj] - rj)
      V[u] <- N * (cov(pji, prj)/nj)
    }
    if (i == s && j != r && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      xr <- samples[[r]]
      ni <- length(xi)
      nj <- length(xj)
      nr <- length(xr)
      ri <- rank(xi)
      rj <- rank(xj)
      rr <- rank(xr)
      rij <- rank(c(xi, xj))
      rir <- rank(c(xi, xr))
      pij <- 1/nj * (rij[1:ni] - ri)
      pir <- 1/nr * (rir[1:ni] - ri)
      V[u] <- -N * (cov(pij, pir)/ni)
    }
    if (i != s && j == r && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      xs <- samples[[s]]
      ni <- length(xi)
      nj <- length(xj)
      ns <- length(xs)
      ri <- rank(xi)
      rj <- rank(xj)
      rs <- rank(xs)
      rji <- rank(c(xj, xi))
      rjs <- rank(c(xj, xs))
      pji <- 1/ni * (rji[1:nj] - rj)
      pjs <- 1/ns * (rjs[1:nj] - rj)
      V[u] <- -N * (cov(pji, pjs)/nj)
    }
  }
  V1 <- matrix(V, ncol = a^2, nrow = a^2)
  switch(effect, weighted = {
    W <- kronecker(t(n/N), diag(a))
    samplesR <- split(rank(response), factorx)
    varsF <- sapply(samplesR,var)
    VH0F = diag(varsF/(n*N^2))
    Si2 = unlist(lapply(1:a,function(arg) var(samplesR[[arg]]-intRanks[[arg]])))
    dfBF = (sum(Si2/(N-n)))^2 / sum((Si2/(N-n))^2/(n-1)) 
    text.output.W <- paste("Global Ranks")
    varKW=sum(c(unlist(samplesR)-(N+1)/2)^2)/(N-1) 
  }, 
  unweighted = {
    W <- kronecker(t(rep(1/a, a)), diag(a))
    samplesR = lapply(1:a,function(arg){
      helpmat=rbind(1:a,matrix(1:a,nrow=a,ncol=a))
      x1 <- samples[[helpmat[1,arg]]]
      help=0
      for(j in 1:a){
        x2 <- samples[[helpmat[j+1,arg]]]
        help=help+1/length(x2)*(rank(c(x1,x2))[1:length(x1)] - rank(x1))
      }
      N/a*help+1/2})

    varsF <- sapply(samplesR,var)
    VH0F = diag(varsF/(n*N^2))
     Si2 = unlist(lapply(1:a,function(arg) var(samplesR[[arg]]-intRanks[[arg]])))
    dfBF = (sum(Si2/(N-n)))^2 / sum((Si2/(N-n))^2/(n-1))
    text.output.W <- paste("Global Pseudo Ranks")
   varKW=sum(c(unlist(samplesR)-(N+1)/2)^2)/(N-1) 
    
  })
  
  pd = W%*%p
  VV=W%*%V1%*%t(W)

  
  result <- list(pd=pd, VH0F=VH0F,VBF=VV, N=N, n=n,dfATS=dfBF,varKW=varKW,placements=placements)
  return(result)
}

###################################################################
# Function for Confidence Limits
###################################################################

Limits <- function(p,V,alpha,N,CI.method){
  
  switch(CI.method, normal={
  Lower <- p - qnorm(1-alpha/2)/sqrt(N)*sqrt(c(diag(V)))
  Upper <- p + qnorm(1-alpha/2)/sqrt(N)*sqrt(c(diag(V)))},
  logit={Psi <- diag(1/(p*(1-p)))
  VLogit <- Psi%*%V%*%t(Psi)
  Lower <- expit(logit(p)- qnorm(1-alpha/2)/sqrt(N)*sqrt(c(diag(VLogit))))
  Upper <- expit(logit(p)+ qnorm(1-alpha/2)/sqrt(N)*sqrt(c(diag(VLogit))))}
  )
  
  res=cbind(Lower=Lower,Upper=Upper)
  return(res)
}

#################################################################
#Wald Type Statistics
#################################################################
Wald <- function(M,H,V){
  WTS = t(H%*%M)%*%ginv(H%*%V%*%t(H))%*%H%*%M
  dfWTS = rankH(H%*%V%*%t(H))
  pv.WTS = 1-pchisq(WTS,dfWTS)
  res.WTS = c(WTS,dfWTS,pv.WTS)
  return(res.WTS)
}

#################################################################
#ANOVA Type Statistics
#################################################################

ANOVATYP <- function(M,H,V,n){
  
  C <- t(H)%*%ginv(H%*%t(H))%*%H
  spur <- sum(diag(C%*%V))
  D <- diag(C)*diag(ncol(C))
  Lambda <- diag(1/(n-1))
  ATS <- 1/spur*t(M)%*%C%*%M
  df_ATS1 <- spur^2/sum(diag(C%*%V%*%C%*%V))
  df_ATS2 <- spur^2/sum(diag(D%*%D%*%V%*%V%*%Lambda))
  
  pv.ATS <- 1-pf(ATS, df_ATS1, df_ATS2)
  res <- c(ATS, df_ATS1, df_ATS2, pv.ATS)
  return(res)
}

ANOVATYPH0P <- function(M,H,V,n,df){
  C <- t(H)%*%ginv(H%*%t(H))%*%H
  spur <- sum(diag(C%*%V))
  ATS <- sum(n)/spur*t(M)%*%C%*%M
  df_ATS1 <- spur^2/sum(diag(C%*%V%*%C%*%V))
  pv.ATS <- 1-pf(ATS, df_ATS1, df)
  res <- c(ATS, df_ATS1, df, pv.ATS)
  return(res)
}

rankH <-function(A){sum(c(diag(ginv(A)%*%A)))}

KWTEST <- function(M,V,n){
N <- sum(n)
KW <- 1/V*sum(n*(N*M+1/2-(N+1)/2)^2)
pv.KW <- 1-pchisq(KW,length(M)-1)
res <- c(KW, length(M)-1, pv.KW)
 return(res)
}   

MCTP <- function(M,H,V,C,n,pla,alpha, sci.method){
#Compute main effects
N <- sum(n) 
pd.main <- H%*%M
V.main <- H%*%V%*%t(H) 
a.main <- length(pd.main)
 CC <- C
CC.H <- CC%*%H
nc <- nrow(CC.H)
CH.pd <- CC%*%pd.main
CH.V <- CC%*%V.main%*%t(CC)
 CH.R <- cov2cor(CH.V)
 
#-----Compute the degrees of freedom of multivariate t-Approximation-----------#
placements <- pla
a <- length(M)
 tmp1 <- sort(rep(1:a, a))
  tmp2 <- rep(1:a, a)
       dfs <- c()
        for (hhh in 1:nc) {
            cc <- CC.H[hhh, ]
            Yk <- list()
            for (l in 1:a) {
                Yk[[l]] <- 0
                for (i in 1:(a^2)) {
                  r <- tmp1[i]
                  s <- tmp2[i]
                  if (s == l && r != l) {
                    Ykhelp <- placements[[i]]
                    Yk[[l]] <- Yk[[l]] + Ykhelp
                  }
                }
            }
            Ykstern <- list()
            for (l in 1:a) {
                Ykstern[[l]] <- cc[l] * Yk[[l]]
                for (i in 1:(a^2)) {
                  r <- tmp1[i]
                  s <- tmp2[i]
                  if (s == l && r != l) {
                    Yksternhelp <- -cc[r] * placements[[i]]
                    Ykstern[[l]] <- Ykstern[[l]] + Yksternhelp
                  }
                }
            }
            variances <- sapply(Ykstern, var)/(a * n)
            varii2 <- (variances == 0)
            variances[varii2] <- 1/N
            dfs[hhh] <- (sum(variances))^2/sum(variances^2/(n - 
                1))
        }
       

    dfT <- round(max(4, min(dfs)))
 
 #------------------Fisher-Transformation of Effects---------------------------# 
 switch(sci.method, multi.t=  {
 
 
SE.CH.pd <- sqrt(c(diag(CH.V)))
T.vec <- sqrt(N)*CH.pd/SE.CH.pd
 pv.main <-sapply(1:nc,function(arg){ 1-pmvt(-abs(T.vec[arg]), abs(T.vec[arg]),delta=rep(0,nc),corr=CH.R,df=dfT)[1]})
 crit <- qmvt(1-alpha, corr = CH.R, tail = "both", df = dfT)$quantile
 Lower <- CH.pd  - crit/sqrt(sum(n))*SE.CH.pd
 Upper <- CH.pd  + crit/sqrt(sum(n))*SE.CH.pd
 res <- data.frame(Effect=CH.pd, Std.Error = sqrt(c(diag(CH.V))/N), T=T.vec, Lower=Lower, Upper=Upper, p.value = pv.main)
 } ,
 
 fisher={
  Cfisher <- 1/2 * log((1 + CH.pd)/(1 - CH.pd))
  Vfisherdev <- diag(c(1/(1 - CH.pd^2)))
  Vfisher <- Vfisherdev %*% CH.V %*% t(Vfisherdev)
  T.vec.Fisher <- sqrt(N) * Cfisher/sqrt(c(diag(Vfisher)))
   crit <- qmvt(1-alpha, corr = CH.R, tail = "both", df = dfT)$quantile
  pv.main.Fisher <-sapply(1:nc,function(arg){ 1-pmvt(-abs(T.vec.Fisher[arg]), abs(T.vec.Fisher[arg]),delta=rep(0,nc),corr=CH.R,df=dfT)[1]})
  Lower.Fisher1 <- Cfisher - crit/sqrt(N) * sqrt(c(diag(Vfisher))) 
Upper.Fisher1 <- Cfisher + crit/sqrt(N) * sqrt(c(diag(Vfisher)))
Lower.Fisher <- (exp(2 * Lower.Fisher1) - 1)/(exp(2 * Lower.Fisher1) + 1)
Upper.Fisher <- (exp(2 * Upper.Fisher1) - 1)/(exp(2 * Upper.Fisher1) + 1)
res <- data.frame(Effect=CH.pd, Std.Error = sqrt(c(diag(CH.V))/N), T=T.vec.Fisher, Lower=Lower.Fisher, Upper=Upper.Fisher, p.value = pv.main.Fisher)
 })  
 
 
 
 

mctp.result= list(Contrast.Matrix=CC.H, Local.Results=res, Global.Result = data.frame(T0=max(abs(res[,3])), p.value=min(res[,6])), DF=dfT, Quantile=crit)

return(mctp.result)

            }
            
MCTP.H0F <- function(M,H,V,C,n,alpha,sci.method){
#Compute main effects 
N <- sum(n)
pd.main <- H%*%M
V.main <- H%*%V%*%t(H) 
a.main <- length(pd.main)
 CC <- C
CC.H <- CC%*%H
nc.main <- nrow(CC.H)
CH.pd <- CC%*%pd.main
CH.V <- CC%*%V.main%*%t(CC)
CH.R <- cov2cor(CH.V)


#-----Compute the degrees of freedom of multivariate t-Approximation-----------#
a <- length(M)
n.V <- n*c(diag(V))
CC.H2 <- CC.H^2
CC.H4 <- CC.H^4
dfTs = sapply(1:nrow(CC.H),function(arg){
(sum(n.V/n*CC.H2[arg,]))^2/sum(n.V^2*CC.H4[arg,]/(n^2*(n-1)))
})

dfT <- round(max(4, min(dfTs)))

switch(sci.method, multi.t={

SE.CH.pd <- sqrt(c(diag(CH.V)))
T.vec <- CH.pd/SE.CH.pd
pv.main <-sapply(1:nc.main,function(arg){ 1-pmvt(-abs(T.vec[arg]), abs(T.vec[arg]),delta=rep(0,nc.main),corr=CH.R,df=dfT)[1]})
crit <- qmvt(1-alpha, corr = CH.R, tail = "both", df = dfT)$quantile

res <- data.frame(Effect=CH.pd, Std.Error = sqrt(c(diag(CH.V))), T=T.vec, p.value = pv.main)
},

fisher={

Cfisher <- 1/2 * log((1 + CH.pd)/(1 - CH.pd))
 Vfisherdev <- diag(c(1/(1 - CH.pd^2)))
 Vfisher <- Vfisherdev %*% CH.V %*% t(Vfisherdev)
 T.vec.Fisher <- Cfisher/sqrt(c(diag(Vfisher)))
 pv.main.Fisher <-sapply(1:nc.main,function(arg){ 1-pmvt(-abs(T.vec.Fisher[arg]), abs(T.vec.Fisher[arg]),delta=rep(0,nc.main),corr=CH.R,df=dfT)[1]})
 crit <- qmvt(1-alpha, corr = CH.R, tail = "both", df = dfT)$quantile
 res <- data.frame(Effect=CH.pd, Std.Error = sqrt(c(diag(CH.V))), T=T.vec.Fisher, p.value = pv.main.Fisher)
}

)
mctp.result= list(Contrast.Matrix=CC.H, Local.Results=res, Global.Result = data.frame(T0=max(abs(res[,3])), p.value=min(res[,4])), DF=dfT, Quantile=crit)

return(mctp.result)

            }
            
            
            
 
 
           
           




#################################################################
#Logit Transformation
#################################################################

logit <- function(p){
  return(log(p/(1-p)))}
expit<-function(p){return(exp(p)/(1+exp(p)))}



################################################################
# for twosample problems
################################################################




BMstat = function(x,y,nx,ny,method){

  Nxy = nx + ny
  xy = c(x,y)
  rxy = rank(xy)
  rx = rank(x)
  ry= rank(y)
  plx = 1/ny*(rxy[1:nx]-rx)
  ply = 1/nx*(rxy[(nx+1):(Nxy)] -ry)
  vx <- var(plx)
  vy <- var(ply)
  vxy = Nxy*(vx/nx + vy/ny)
  pd=mean(ply)
  pd0=(pd==0)
  pd1 = (pd==1)
  pd[pd0] = (1/(2*nx))/ny
  pd[pd1] = (nx-1/(2*ny))/nx
  vxy0= (vxy==0)
  vxy[vxy0] = Nxy/(2*nx*ny)
  VH0F <- 1/(Nxy*(Nxy-1))*sum((rxy-(Nxy+1)/2)^2)*(1/nx+1/ny)
      df.sw <- (vx/nx + vy/ny)^2/(vx^2/(nx^2*(nx - 1)) + vy^2/(ny^2*(ny - 1)))
  df.sw[is.nan(df.sw)] <- 1000
  Z.exact=sum(rxy[(nx+1):Nxy])
  

  switch(method,
  normal={
  T = sqrt(Nxy)*(pd-1/2)/sqrt(vxy)
  se.pd=sqrt(vxy/Nxy)
  f.pd=pd
  },
  logit ={
  slogit=(1/(pd*(1-pd))*sqrt(vxy))
  T = sqrt(Nxy)*logit(pd)/(slogit)
  se.pd=slogit/sqrt(Nxy)
  f.pd=logit(pd)
  },
  probit={
  sprobit=sqrt(2*pi)/(exp(-1/2*(qnorm(pd))^2))*sqrt(vxy)
  T= sqrt(Nxy) * qnorm(pd)/(sprobit)
  se.pd=sprobit/sqrt(Nxy)
  f.pd=qnorm(pd)
  },
  t.app={
  T = sqrt(Nxy)*(pd-1/2)/sqrt(vxy)
  se.pd=sqrt(vxy/Nxy)
  f.pd=pd
  }
  )
  
  res <-data.frame(pd=pd, vx=vx, vy=vy, vxy=vxy, T=T, f.pd =f.pd, se.pd=se.pd,VH0F=VH0F,df.sw=df.sw,se.output=sqrt(vxy/Nxy),Z.exact=Z.exact)
  return(res)
}

rev.function=function(x,method){
switch(method, normal={identity(x)},
t.app={identity(x)},
logit={expit(x)},
probit=pnorm(x))
}


limits2 <-function(f.pd,crit1,crit2,se.pd,alternative,method){

switch(alternative,
two.sided={
Lower=rev.function(f.pd-crit1*se.pd,method)
Upper=rev.function(f.pd-crit2*se.pd,method)},

less={
Lower= 0
Upper = rev.function(f.pd+crit2*se.pd,method)
},

greater={
Lower = rev.function(f.pd-crit2*se.pd,method)
Upper=1
}

)
res =data.frame(Lower=Lower,Upper=Upper)
return(res)
}




