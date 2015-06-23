
# I generate some data here, but the process starts with any
# species trinary matrix
require(binom)


#Function that creates some fake data for a species
#that is rare
f.data.generator<-function(sites,days,psi,p,phi,gamma,nyears) {
  #first year of data
  y1<-matrix(NA,nr=sites,nc=days)
  #generate the expected occupancies
  z<-rbinom(sites,1,psi)
  #generate the observations
  for(i in 1:sites)
    y1[i,]<-rbinom(days,1,z[i]*p)
  #subsequent years
  #three dimensional matrix to store the results
  yk<-array(NA,dim=c(sites,days,nyears))
  yk[,,1]<-y1
  for(k in 2:nyears){
    #generate the deterministic part of the model
    occ<-apply(yk[,,k-1],1,max,na.rm=T)
    z<-rbinom(sites,1,occ*phi+(1-occ)*gamma)
    #generate the observations
    for(i in 1:sites)
      yk[i,,k]<-rbinom(days,1,z[i]*p)
    
  }  
  yk
}

#function to calculate the mode of a distribution
f.mode<-function(data){
  qwe<-density(data)
  qwe$x[which(qwe$y==max(qwe$y))]
  
}

# Create a new data set
 


data<-f.data.generator(sites=60,days=15,psi=0.05,p=.05,phi=0.4,gamma=0.1,nyears=5)
#simulating some holes in the data
data[31:60,1:7,]<-NA
data[1:30,8:15,]<-NA
data[60,,]<-NA
data[13,,]<-NA

nyears <- dim(data)[3]
# calculate whether the species was present or absent at each point each year
pa.year<-apply(data,c(1,3),max,na.rm=T)
pa.year[pa.year=="-Inf"]<-NA

# how many camera traps were functional?
n.cams <- apply(pa.year,2,function(x) {sum(!is.na(x))})

# How many detections of each species per year?
n.dets <- apply(pa.year,2,sum,na.rm=T)

# Average number of detections? If mean.dets < 5 we use this process to calculate
# occupancy
mean.dets <- mean(n.dets)

# naive occupancy
occ <- n.dets/n.cams
#new.occ will contain 1000 realizations of each 
new.occ <- matrix(NA,nr=1000,nc=nyears)

for(i in 1:nyears){
  conf.bin <- binom.bayes(x=n.dets[i],n=n.cams[i],conf.level=0.95,type="highest")
  mod.bin <- binom.bayes.densityplot(bayes=conf.bin)
  new.occ[,i] <- sample(x=mod.bin$data$xx,size=1000,replace=T,prob=mod.bin$data$yy)
}

# Check everything makes sense and graph it
conf.lims <- apply(new.occ,2,quantile,c(0.025,0.975))
occ.mode <- apply(new.occ,2,f.mode)

#plot mode of the distribution plus confidence limits
plot(1:nyears,occ.mode,t='l',ylim=c(0,1))
lines(1:nyears,conf.lims[1,],lty=2)
lines(1:nyears,conf.lims[2,],lty=2)
#observed
points(1:nyears,occ)

