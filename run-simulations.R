#Execute simulations with different parameter values
require(R2jags)
require(ggplot2)
require(reshape2)
system.time(res<-f.fit.JAGS(psi=0.01,p=0.01,phi=0.01,gamma=0.01,nsites=60,ndays=15,nyears=5,ns=1))
#run one with 500 sumulations
system.time(res.01.01.500<-f.fit.JAGS(psi=0.01,p=0.01,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=500))
system.time(res.01.05<-f.fit.JAGS(psi=0.01,p=0.05,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=100))
system.time(res.01.1<-f.fit.JAGS(psi=0.01,p=0.1,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=100))
system.time(res.01.5<-f.fit.JAGS(psi=0.01,p=0.5,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=100))

#psi=0.5, p=0.01,0.05,0.1,0.5, 0.8
system.time(res.5.01<-f.fit.JAGS(psi=0.5,p=0.01,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=100))
save(res.5.01,file="res.5.01")
system.time(res.5.05<-f.fit.JAGS(psi=0.5,p=0.05,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=100))
save(res.5.05,file="res.5.05")
system.time(res.5.1<-f.fit.JAGS(psi=0.5,p=0.1,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=100))
save(res.5.1,file="res.5.1")
system.time(res.5.5<-f.fit.JAGS(psi=0.5,p=0.5,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=100))
save(res.5.5,file="res.5.5")
system.time(res.5.8<-f.fit.JAGS(psi=0.5,p=0.8,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=3,ns=100))
save(res.5.8,file="res.5.8")

#check if results fit the expected patterns
qwe<-melt(res.01.01)
qwe$L1<-as.factor(qwe$L1)
names(qwe)<-c("parm","value","stat")
psi.01.01<-qwe[qwe$parm=="psi",]
qplot(stat,value,data=psi.01.01,geom="boxplot")
p<-ggplot(psi.01.01,aes(x=stat,y=value))
p<-p+geom_boxplot()+geom_abline(inter=0.01,slope=0)

#displaying all variables
res.01.01<-qwe
p<-ggplot(res.01.01,aes(x=stat,y=value))
p<-p+geom_boxplot()+geom_abline(inter=0.1,slope=0,subset=.(parm=="psi"))
p

#using faceting to display varables
p<-ggplot(res.01.01,aes(x=stat,y=value))
p<-p+facet_grid(parm~.,scales="free_y")+geom_jitter(alpha=0.5,size=1.5)+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_abline(inter=0.01,slope=0,data=subset(res.01.01,parm=="psi"))+geom_abline(inter=0.01,slope=0,data=subset(res.01.01,parm=="p"))+geom_abline(inter=0.1,slope=0,data=subset(res.01.01,parm=="phi" | parm=="gamma"))
p
ggsave(file="res.01.01.pdf",width=10,height=18,units="cm",plot=p)
#check the other one
qwe<-melt(res.01.05)
qwe$L1<-as.factor(qwe$L1)
names(qwe)<-c("parm","value","stat")
res.01.05<-qwe

p<-ggplot(res.01.05,aes(x=stat,y=value))
p<-p+facet_grid(parm~.,scales="free_y")+geom_boxplot()+geom_abline(inter=0.01,slope=0,subset=.(parm=="psi"))+geom_abline(inter=0.05,slope=0,subset=.(parm=="p"))+geom_abline(inter=0.1,slope=0,subset=.(parm=="phi" | parm=="gamma"))
p

#check the run with 500 sinulations
qwe<-melt(res.01.01.500)
qwe$L1<-as.factor(qwe$L1)
names(qwe)<-c("parm","value","stats")
res.01.01.500<-qwe

p<-ggplot(res.01.01.500,aes(x=stats,y=value))
p<-p+facet_grid(parm~.,scales="free_y")+geom_jitter(alpha=0.5,size=1.5)+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_abline(inter=0.01,slope=0,data=subset(res.01.01.500,parm=="psi"))+geom_abline(inter=0.01,slope=0,data=subset(res.01.01.500,parm=="p"))+geom_abline(inter=0.1,slope=0,data=subset(res.01.01.500,parm=="phi" | parm=="gamma"))
p
ggsave(file="res.01.01.500.pdf",width=10,height=18,units="cm",plot=p)


#After simulations were finished now put the whole thing in one large data frame for ggplot

#the first one had gone through some transformations so easy to put in the format
res.01.01<-data.frame(res.01.01,psi=0.01,p=0.01)
#for all the other ones
#first load the simulation from file
load("Simulation results/res.01.05")
#melt
res.01.05<-melt(res.01.05)
#add psi and p as variables
res.01.05<-data.frame(res.01.05,psi=0.01,p=0.05)
#change variable names around
res.01.05$L1<-as.factor(res.01.05$L1)
names(res.01.05)[1:3]<-c("parm","value","stats")

#do a little function to do this

#f.ConvertListToDataframe<-function(list,psi,p){
res.01.01<-f.ConvertListToDataframe(res.01.01,.01,.01)
res.01.05<-f.ConvertListToDataframe(res.01.05,.01,.05)
res.01.1<-f.ConvertListToDataframe(res.01.1,.01,.1)
#res.01.5<-f.ConvertListToDataframe(res.01.5,.01,.5)
#res.01.8<-f.ConvertListToDataframe(res.01.8,.01,.8)

#put all the ones with the same psi together
res.01<-rbind(res.01.01,res.01.05,res.01.1)#,res.01.5,res.01.8)
#for psi=0.1
res.1.01<-f.ConvertListToDataframe(res.1.01,.1,.01)
res.1.05<-f.ConvertListToDataframe(res.1.05,.1,.05)
res.1.1<-f.ConvertListToDataframe(res.1.1,.1,.1)
#res.1.5<-f.ConvertListToDataframe(res.1.5,.1,.5)
#res.1.8<-f.ConvertListToDataframe(res.1.8,.1,.8)

res.1<-rbind(res.1.01,res.1.05,res.1.1)#,res.1.5,res.1.8)

#for psi=0.05
res.05.01<-f.ConvertListToDataframe(res.05.01,.05,.01)
res.05.05<-f.ConvertListToDataframe(res.05.05,.05,.05)
res.05.1<-f.ConvertListToDataframe(res.05.1,.05,.1)
#res.05.5<-f.ConvertListToDataframe(res.05.5,.05,.5)
#res.05.8<-f.ConvertListToDataframe(res.05.8,.05,.8)

res.05<-rbind(res.05.01,res.05.05,res.05.1)#,res.05.5,res.05.8)

#for psi=0.5
res.5.01<-f.ConvertListToDataframe(res.5.01new,.5,.01)
res.5.05<-f.ConvertListToDataframe(res.5.05,.5,.05)
res.5.1<-f.ConvertListToDataframe(res.5.1,.5,.1)
#res.5.5<-f.ConvertListToDataframe(res.5.5,.5,.5)
#res.5.8<-f.ConvertListToDataframe(res.5.8,.5,.8)

res.5<-rbind(res.5.01,res.5.05,res.5.1)#,res.5.5,res.5.8)


#for psi=0.8
res.8.01<-f.ConvertListToDataframe(res.8.01,.8,.01)
res.8.05<-f.ConvertListToDataframe(res.8.05,.8,.05)
res.8.1<-f.ConvertListToDataframe(res.8.1,.8,.1)
#res.8.5<-f.ConvertListToDataframe(res.8.5,.8,.5)
#res.8.8<-f.ConvertListToDataframe(res.8.8,.8,.8)

res.8<-rbind(res.8.01,res.8.05,res.8.1)#,res.8.5,res.8.8)

#put everything together
res<-rbind(res.01,res.05,res.1,res.5,res.8)


# do a graph for psi
res.psi<-subset(res,parm=="psi")
res.psi<-droplevels(res.psi)
p<-ggplot(data=res.psi,aes(x=stats,y=value))

p<-p+facet_grid(psi~p,scales="free_y")+geom_violin(fill="blue",alpha=0.1)+geom_boxplot(width=0.1,outlier.size=0,fill='grey50')+stat_summary(fun.y=median,geom="point",fill="white",size=4,shape=21)
p<-p+geom_abline(inter=0.01,slope=0,data=subset(res.psi,psi==0.01))
p<-p+geom_abline(inter=0.05,slope=0,data=subset(res.psi,psi==0.05))
p<-p+geom_abline(inter=0.1,slope=0,data=subset(res.psi,psi==0.1))
p<-p+geom_abline(inter=0.5,slope=0,data=subset(res.psi,psi==0.5))
p<-p+geom_abline(inter=0.8,slope=0,data=subset(res.psi,psi==0.8))
p<-p+labs(x="Measure of central tendecy",y="Recovered initial occupancy")
gpsi<-p
pdf("test-psi2.pdf",width=15/2.54,height=20/2.54)
#pushViewport(viewport(layout=grid.layout(2,2)))
print(gpsi,vp=viewport(width=0.9,height=0.9))
grid.text(0.95,0.5,label = "Model initial occupancy", rot = 270)
grid.text(0.5,0.95,label= "Model detection probability")
#print(gp,vp=viewport(layout.pos.row=2,layout.pos.col=1))
#print(ggamma,vp=viewport(layout.pos.row=1,layout.pos.col=2))
#print(gphi,vp=viewport(layout.pos.row=2,layout.pos.col=2))
dev.off()

#ggsave(file="res.psi.pdf",width=15,height=20,units="cm",plot=p)

# do a graph for p

p<-ggplot(data=subset(res,parm=="p"),aes(x=stats,y=value))

p<-p+facet_grid(p~psi,scales="free_y")+geom_violin(fill="blue",alpha=0.1)+geom_boxplot(width=0.1,outlier.size=0,fill='grey50')+stat_summary(fun.y=median,geom="point",fill="white",size=4,shape=21)
#+geom_jitter(alpha=0.5,size=1.5)
p<-p+geom_abline(inter=0.01,slope=0)
p<-p+geom_abline(inter=0.05,slope=0,data=subset(res,p==.05))
p<-p+geom_abline(inter=0.1,slope=0,data=subset(res,p==.1))
p<-p+labs(x="Measure of central tendecy",y="Recovered detection estimate",title="Model Occupancy",size=2)
gp<-p
ggsave(file="res.p.pdf",width=15,height=20,units="cm",plot=p)

#p<-p+geom_abline(inter=0.5,slope=0,data=subset(res,p==.5))
#p<-p+geom_abline(inter=0.8,slope=0,data=subset(res,p==.8))
#p<-p+labs(x="Measure of central tendecy",y="p ",title="Detection probability estimate")

#+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_abline(inter=0.01,slope=0,data=subset(res.01.01.500,parm=="psi"))+geom_abline(inter=0.01,slope=0,data=subset(res.01.01.500,parm=="p"))+geom_abline(inter=0.1,slope=0,data=subset(res.01.01.500,parm=="phi" | parm=="gamma"))
p
ggsave(file="res.p.pdf",width=25,height=15,units="cm",plot=p)

#do a graph for gamma
p<-ggplot(data=subset(res,parm=="gamma"),aes(x=stats,y=value))
p<-p+facet_grid(psi~p,scales="free_y")+geom_violin(fill="blue",alpha=0.1)+geom_boxplot(width=0.1,outlier.size=0,fill='grey50')+stat_summary(fun.y=median,geom="point",fill="white",size=4,shape=21)
p<-p+geom_hline(yintercept=0.1)
p<-p+labs(x="Measure of central tendecy",y="Recovered estimate of gamma",title="Model detection probability")
ggamma<-p
ggsave(file="res.gamma.pdf",width=15,height=20,units="cm",plot=ggamma)

#do a graph for phi
p<-ggplot(data=subset(res,parm=="phi"),aes(x=stats,y=value))
p<-p+facet_grid(psi~p,scales="free_y")+geom_violin(fill="blue",alpha=0.1)+geom_boxplot(width=0.1,outlier.size=0,fill='grey50')+stat_summary(fun.y=median,geom="point",fill="white",size=4,shape=21)
p<-p+geom_hline(yintercept=0.2)
p<-p+labs(x="Measure of central tendecy",y="Recovered estimate of phi",title="Model detection probability")
gphi<-p
ggsave(file="res.phi.pdf",width=15,height=20,units="cm",plot=gphi)


#check the run of psi=0.05, p=0.01 with 500 runs
p<-ggplot(res.05.01.500,aes(x=stats,y=value))
p<-p+facet_grid(parm~.,scales="free_y")+geom_jitter(alpha=0.5,size=1.5)+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_abline(inter=0.05,slope=0,data=subset(res.05.01.500,parm=="psi"))+geom_abline(inter=0.01,slope=0,data=subset(res.05.01.500,parm=="p"))+geom_abline(inter=0.1,slope=0,data=subset(res.05.01.500,parm=="phi" | parm=="gamma"))
p
ggsave(file="res.05.01.500.pdf",width=10,height=18,units="cm",plot=p)

#check the run of psi=0.05, p=0.01 with 100 runs and 30 days (instead of 15)
p<-ggplot(res.05.01.30d,aes(x=stats,y=value))
p<-p+facet_grid(parm~.,scales="free_y")+geom_jitter(alpha=0.5,size=1.5)+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_abline(inter=0.05,slope=0,data=subset(res.05.01.30d,parm=="psi"))+geom_abline(inter=0.01,slope=0,data=subset(res.05.01.30d,parm=="p"))+geom_abline(inter=0.1,slope=0,data=subset(res.05.01.30d,parm=="phi" | parm=="gamma"))
p
ggsave(file="res.05.01.30s-noPriors.pdf",width=10,height=18,units="cm",plot=p)

#check the run of psi=0.05, p=0.01 with 100 runs and 15 days (no priors)
p<-ggplot(res.01.01.15d,aes(x=stats,y=value))
p<-p+facet_grid(parm~.,scales="free_y")+geom_jitter(alpha=0.5,size=1.5)+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_abline(inter=0.01,slope=0,data=subset(res.01.01.15d,parm=="psi"))+geom_abline(inter=0.01,slope=0,data=subset(res.01.01.15d,parm=="p"))+geom_abline(inter=0.1,slope=0,data=subset(res.01.01.15d,parm=="phi" | parm=="gamma"))
p
ggsave(file="res.05.01.15d-noPriors.pdf",width=10,height=18,units="cm",plot=p)

#look at the effect of number of camera traps
#90 cameras
res.01.01n90<-f.ConvertListToDataframe(res.01.01n90,.01,.01)
p<-ggplot(res.01.01n90,aes(x=stats,y=value))
p<-p+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_hline(yintercept=0.01)

#30 cameras
res.01.01n30<-f.ConvertListToDataframe(res.01.01n30,.01,.01)
p<-ggplot(res.01.01n30,aes(x=stats,y=value))
p<-p+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_hline(yintercept=0.01)
p

#60 points but 100 days
res.5.01new<-f.ConvertListToDataframe(res.5.01new,.5,.01)
p<-ggplot(res.5.01new,aes(x=stats,y=value))
p<-p+geom_boxplot(fill="blue",alpha=0.1,outlier.size=1)+geom_hline(yintercept=0.5)
p


#running it with 100 days
res.5.01new<-f.ConvertListToDataframe(res.5.01new,.5,.01)
p<-ggplot(res.5.01new,aes(x=stats,y=value))
p<-p+facet_grid(parm~.,scales="free_y")+geom_violin(fill="blue",alpha=0.1)+geom_boxplot(width=0.1,outlier.size=0,fill='grey50')+stat_summary(fun.y=median,geom="point",fill="white",size=4,shape=21)+geom_abline(inter=0.5,slope=0,data=subset(res.5.01new,parm=="psi"))+geom_abline(inter=0.01,slope=0,data=subset(res.5.01new,parm=="p"))+geom_abline(inter=0.1,slope=0,data=subset(res.5.01new,parm=="phi" | parm=="gamma"))
p
ggsave(file="res.5.01-100days.pdf",width=10,height=18,units="cm",plot=p)

#running it with 50 days instead of 100
res.5.01new2<-f.ConvertListToDataframe(res.5.01new2,.5,.01)
p<-ggplot(res.5.01new2,aes(x=stats,y=value))
p<-p+facet_grid(parm~.,scales="free_y")+geom_violin(fill="blue",alpha=0.1)+geom_boxplot(width=0.1,outlier.size=0,fill='grey50')+stat_summary(fun.y=median,geom="point",fill="white",size=4,shape=21)+geom_abline(inter=0.5,slope=0,data=subset(res.5.01new2,parm=="psi"))+geom_abline(inter=0.01,slope=0,data=subset(res.5.01new2,parm=="p"))+geom_abline(inter=0.1,slope=0,data=subset(res.5.01new2,parm=="phi"))+geom_abline(inter=0.15,slope=0,data=subset(res.5.01new2,parm=="gamma"))
p

ggsave(file="res.5.01-50days.pdf",width=10,height=18,units="cm",plot=p)

require(grid)
pdf("All parameters.pdf",width=25,height=25)
pushViewport(viewport(layout=grid.layout(2,2)))
print(gpsi,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(gp,vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(ggamma,vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(gphi,vp=viewport(layout.pos.row=2,layout.pos.col=2))
dev.off()

#run a representative realization and graph it
realz<-f.data.generator(sites=60,days=100,psi=0.01,p=0.01,phi=0.2,gamma=0.1,nyears=3)
jags.data <- list(y = realz, nsite = 60, nrep = 100, nyear = 3)

# Initial values
initial <- apply(realz, c(1,3), max, na.rm = TRUE)
#tmp[tmp=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
inits <- function(){ list(z = initial)}
#params to monitor
params <- c("psi", "phi", "gamma", "p", "lambda") 
# MCMC settings
ni <- 8000
nt <- 2
nb <- 3000
nc <- 4

out <- jags(jags.data, inits, params, "Dynocc-jags-psiYear1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Epsi<-out$BUGSoutput$sims.list$psi
Mpsi<-f.cal.psis(0.01,0.2,0.1,3)
p<-graph.psi(psi=Epsi,Mpsi,fun=f.mode,"")
ggsave(file="Figure5.pdf",width=10,height=6,units="cm",plot=p)

# do it with psi=0.05
realz<-f.data.generator(sites=60,days=100,psi=0.05,p=0.01,phi=0.2,gamma=0.1,nyears=3)

initial <- apply(realz, c(1,3), max, na.rm = TRUE)
#tmp[tmp=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
inits <- function(){ list(z = initial)}
jags.data <- list(y = realz, nsite = 60, nrep = 100, nyear = 3)
out <- jags(jags.data, inits, params, "Dynocc-jags-psiYear1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Epsi2<-out$BUGSoutput$sims.list$psi
Mpsi2<-f.cal.psis(0.05,0.2,0.1,3)
p<-graph.psi(psi=Epsi2,Mpsi2,fun=f.mode,"")

# do it with psi=0.1
realz<-f.data.generator(sites=60,days=100,psi=0.1,p=0.01,phi=0.2,gamma=0.1,nyears=3)

initial <- apply(realz, c(1,3), max, na.rm = TRUE)
#tmp[tmp=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
inits <- function(){ list(z = initial)}
jags.data <- list(y = realz, nsite = 60, nrep = 100, nyear = 3)
out <- jags(jags.data, inits, params, "Dynocc-jags-psiYear1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Epsi3<-out$BUGSoutput$sims.list$psi
Mpsi3<-f.cal.psis(0.1,0.2,0.1,3)
p<-graph.psi(psi=Epsi3,Mpsi3,fun=f.mode,"")