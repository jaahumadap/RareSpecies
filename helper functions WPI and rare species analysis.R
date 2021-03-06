#HELPER functions in WPI analysis and rare species analyses
#create some fake data for a species
#that is rare
mean.det.prob <- function(spmat){
  #calculate the observed mean det probability for a species
  mean(apply(spMat,c(3),sum,na.rm=T)/apply(spMat,c(3),function(x){sum(!is.na(x))}))
}

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
#function to recursively calculate psi for year 2,3, etc.

f.cal.psis<-function(psi1,phi,gamma,nyears){
	res<-numeric()
	res[1]<-psi1
	for(i in 2:nyears)
		res[i]<-res[i-1]*phi+(1-res[i-1])*gamma
	res
}
#function to calculate the mode of a distribution
f.mode<-function(data){
	qwe<-density(data)
	qwe$x[which(qwe$y==max(qwe$y))]
	
}


#Main function to execute an iteration of the fitting process with JAGS

f.fit.JAGS<-function(psi=0.1,p=0.1,phi=0.1,gamma=0.1,nsites=60,ndays=15,nyears=5,nsims=1){
  #set up data frames to store results
  res.mean<-data.frame(psi=rep(NA,nsims),p=NA,phi=NA,gamma=NA)
  res.median<-data.frame(psi=rep(NA,nsims),p=NA,phi=NA,gamma=NA)
  res.mode<-data.frame(psi=rep(NA,nsims),p=NA,phi=NA,gamma=NA)
  res.diags <- data.frame(ndets=rep(NA,nsims),ndetsy=rep(NA,nsims),diff=rep(NA,nsims),fit.score=rep(NA,nsims))
  for(i in 1:nsims){
    print(paste("Simulation ",i," of ",nsims,sep=""))
    #generate some data first
    data<-f.data.generator(nsites,ndays,psi,p,phi,gamma,nyears)
    # n of pres/abs per year
    n.dets <- apply(apply(data,c(1,3),max,na.rm=T),2,sum)
    #n detections per site per year
    n.dets.sy <- mean(apply(apply(data,c(1,3),sum),2,sum)/nsites)
    naive.occ <-n.dets/nsites
    
    #setup data for JAGS
    jags.data <- list(y = data, nsite = nsites, nrep = ndays, nyear = nyears)
    
    # Initial values
    initial <- apply(data, c(1,3), max, na.rm = TRUE)
    inits <- function(){ list(z = initial)}
    #params to monitor
    params <- c("psi", "phi", "gamma", "p", "lambda") 
    # MCMC settings
    ni <- 2000
    nt <- 2
    nb <- 1000
    nc <- 3
    
    out <- jags(jags.data, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
  
    rawSims<-out$BUGSoutput$sims.array
    pnames<-dimnames(rawSims)[[3]]
    rawSims<-rawSims[,,c(2,5:9)]
    dimnames(rawSims)<-list(NULL,NULL,names=pnames[c(2,5:9)])
    
    #calculate mean, median and mode and store in a dataframe
    res.mean$psi[i]<-mean(rawSims[,,"psi[1]"],na.rm=T)
    res.mean$p[i]<-mean(rawSims[,,"p"],na.rm=T)
    res.mean$phi[i]<-mean(rawSims[,,"phi"],na.rm=T)
    res.mean$gamma[i]<-mean(rawSims[,,"gamma"],na.rm=T)
    
    res.median$psi[i]<-median(rawSims[,,"psi[1]"],na.rm=T)
    res.median$p[i]<-median(rawSims[,,"p"],na.rm=T)
    res.median$phi[i]<-median(rawSims[,,"phi"],na.rm=T)
    res.median$gamma[i]<-median(rawSims[,,"gamma"],na.rm=T)
    
    res.mode$psi[i]<-f.mode(rawSims[,,"psi[1]"])
    res.mode$p[i]<-f.mode(rawSims[,,"p"])
    res.mode$phi[i]<-f.mode(rawSims[,,"phi"])
    res.mode$gamma[i]<-f.mode(rawSims[,,"gamma"])
    
    #where are the psi's in the JAGS object
    indx <- grep(pattern="psi",rownames(out$BUGSoutput$summary))
    #model fits
    fitted.psi <- out$BUGSoutput$summary[indx,2]
    low.psi <- out$BUGSoutput$summary[indx,3]
    hi.psi <- out$BUGSoutput$summary[indx,7]
    #mean difference between observed and expected
    res.diags$diff[i] <- sum(abs(fitted.psi-naive.occ))/nyears
    #the fit score has a max value of nyears if observed occupancy is within the boundaries of the model
    res.diags$fit.score[i] <- sum(naive.occ < hi.psi & naive.occ > low.psi)
    res.diags$ndets[i] <-mean(n.dets)
    res.diags$ndetsy[i] <- n.dets.sy
  }
  list(mean=res.mean,median=res.median,mode=res.mode,diags=res.diags)
}

f.ConvertListToDataframe<-function(list,psi,p){
  
  load(paste("New simulations/",deparse(substitute(list)),sep=""),.GlobalEnv)
  #melt
  df<-melt(list)
  #add psi and p as variables
  df<-data.frame(df,psi=psi,p=p)
  #change variable names around
  #df$L1<-as.factor(df$L1)
  names(df)[1:3]<-c("parm","value","stats")
  df
}

#another version of the function without passing priors to jags
f.fit.JAGS2<-function(psi=0.05,p=0.05,phi=0.2,gamma=0.2,nsites=30,ndays=10,nyears=1,nsims=1){
  #set up data frames to store results
  res.mean<-data.frame(psi=rep(NA,nsims),p=NA,phi=NA,gamma=NA)
  res.median<-data.frame(psi=rep(NA,nsims),p=NA,phi=NA,gamma=NA)
  res.mode<-data.frame(psi=rep(NA,nsims),p=NA,phi=NA,gamma=NA)
  for(i in 1:nsims){
    print(paste("Simulation ",i," of ",nsims,sep=""))
    #generate some data first
    data<-f.data.generator(nsites,ndays,psi,p,phi,gamma,nyears)
    
    #estimate detection probability from the data
    #ndetections<-apply(data,c(1,3),sum,na.rm=T)
    #nsamples<-apply(!is.na(data),c(1,3),sum)
    #etProb<-ndetections/nsamples
    #indx<-which(detProb>0)
    #if there are no 1's in the data
    #if(!length(indx)){
    #  muDetPrior<-0.01
    #} else
    #  muDetPrior<-mean(detProb[indx])
    
    #estimate occupancy in year1 from the data
    #muPsi1Prior<-sum(apply(data[,,1],1,max,na.rm=T))/nsites
    #if there are no 1's in the data
    #if(!muPsi1Prior){
    #  muPsi1Prior<-0.01
    #}
    
    #setup data for JAGS
    jags.data <- list(y = data, nsite = nsites, nrep = ndays, nyear = nyears)
    
    # Initial values
    initial <- apply(data, c(1,3), max, na.rm = TRUE)
    #tmp[tmp=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
    inits <- function(){ list(z = initial)}
    #params to monitor
    params <- c("psi", "phi", "gamma", "p", "lambda") 
    # MCMC settings
    ni <- 8000
    nt <- 2
    nb <- 3000
    nc <- 4
    
    out <- jags(jags.data, inits, params, "Dynocc-jags-psiYear1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
    
    rawSims<-out$BUGSoutput$sims.array
    pnames<-dimnames(rawSims)[[3]]
    rawSims<-rawSims[,,c(2,5:9)]
    dimnames(rawSims)<-list(NULL,NULL,names=pnames[c(2,5:9)])
    
    #calculate mean, median and mode and store in a dataframe
    res.mean$psi[i]<-mean(rawSims[,,"psi[1]"],na.rm=T)
    res.mean$p[i]<-mean(rawSims[,,"p"],na.rm=T)
    res.mean$phi[i]<-mean(rawSims[,,"phi"],na.rm=T)
    res.mean$gamma[i]<-mean(rawSims[,,"gamma"],na.rm=T)
    
    res.median$psi[i]<-median(rawSims[,,"psi[1]"],na.rm=T)
    res.median$p[i]<-median(rawSims[,,"p"],na.rm=T)
    res.median$phi[i]<-median(rawSims[,,"phi"],na.rm=T)
    res.median$gamma[i]<-median(rawSims[,,"gamma"],na.rm=T)
    
    res.mode$psi[i]<-f.mode(rawSims[,,"psi[1]"])
    res.mode$p[i]<-f.mode(rawSims[,,"p"])
    res.mode$phi[i]<-f.mode(rawSims[,,"phi"])
    res.mode$gamma[i]<-f.mode(rawSims[,,"gamma"])
  }
  list(mean=res.mean,median=res.median,mode=res.mode)
}
