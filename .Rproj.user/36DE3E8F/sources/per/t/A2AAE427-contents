
psr=function(formula, data, psranks){
 dat.Model0 <- model.frame(formula, data)
dat.Model0[,2] = as.factor(dat.Model0[,2])
dat.Model0 <- dat.Model0[order(dat.Model0[, 2]), ]
n <- aggregate(formula, data = dat.Model0, length)[,2]
gn=sum(n)
a=length(n)
n1     = n[1]                                                            
vsum   = matrix(1/n1,nrow=n1,ncol=1)                                                 
for(i in 2:a){                                                             
   ni = n[i]                                                             
 vsum = rbind(vsum, matrix(1/ni, nrow=ni,ncol=1) )                                    
}

X=dat.Model0[,1]
r=X
gn=length(X)
pr=c()
for ( k in 1:gn){                                                            
   pr[k] = ((sign(r[k]*matrix(1,nrow=1,ncol=gn) - t(r)) + matrix(1,nrow=1,ncol=gn))*0.5)%*%vsum*gn/a + 0.5 
} 

dat.Model0$psr=pr
names(dat.Model0)[3] = psranks
dat.Model0
}

