
# I generate some data here, but the process starts with any
# species trinary matrix
require(binom)
data<-f.data.generator(sites=60,days=15,psi=.01,p=.01,phi=0.05,gamma=0.05,nyears=5)
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
occ <- n.dets/n
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
