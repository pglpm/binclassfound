## Author: PGL  Porta Mana
## Created: 2022-03-17T14:21:57+0100
## Last-Updated: 2022-03-17T15:31:07+0100
################
## Exploration of several issues for binary classifiers
################
if(file.exists("/cluster/home/pglpm/R")){
    .libPaths(c("/cluster/home/pglpm/R",.libPaths()))
}
#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
## library('khroma')
## palette(colour('bright')())
## scale_colour_discrete <- scale_colour_bright
## palette(colour('muted')())
library('data.table')
## library('ggplot2')
## library('ggthemes')
## theme_set(theme_bw(base_size=18))
#library('cowplot')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
print('availableCores:')
print(availableCores())
print('availableCores-multicore:')
print(availableCores('multicore'))
if(file.exists("/cluster/home/pglpm/R")){
    plan(multicore, workers=availableCores()-1)
}else{
    plan(multisession, workers=6)
}
##library('ash')
## library('LaplacesDemon')
## library('extraDistr')
## library('mvtnorm')
## options(bitmapType='cairo')
## pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
## pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
## library('nimble')
#### End custom setup ####

condpRows <- function(x){
    t(t(x)/colSums(x))
}
condpCols <- function(x){
    x/rowSums(x)
}

## We consider two binary {0,1} variables X and Y and a population with a given conditional frequency of Y given X

pXgY0 <- rbind(0.375, 1-0.375)
pXgY1 <- rbind(1-0.25, 0.25)
##
pXgY <- rbind(t(pXgY0), t(pXgY1))
##
pY <- rbind(1-0.125, 0.125)
##
pYX <- pXgY * c(pY)
pX <- t(pXgY) %*% pY
##  0.372, 0.628
sum((pXgY - condpCols(pYX))^2)
pYgX <- condpRows(pYX)
##           [,1]       [,2]
## [1,] 0.7777778 0.94594595
## [2,] 0.2222222 0.05405405


pY2 <- rbind(0.5, 1-0.5)
##
pYX2 <- pXgY * c(pY2)
pX2 <- t(pXgY) %*% pY2
##  0.372, 0.628
sum((pXgY - condpCols(pYX2))^2)
pYgX2 <- condpRows(pYX2)
##           [,1]      [,2]
## [1,] 0.3333333 0.7142857
## [2,] 0.6666667 0.2857143



#######################################
#### FUNCTION TO CALCULATE MUTUAL INFO FROM JOINT DISTRIBUTION
## freqs[S,B] = freq spike count B and stimulus S (one ROW per stimulus)
## The function calculates the conditional frequencies of B|S
## and if requested constructs
## a new joint distribution with equal marginals for S
## Note: don't need to normalize input to mutualinfo
myoptim <- function(par, fn){
    resu0 <- list(par=par)
    resu <- optim(par=resu0$par, fn=fn, control=list(factr = 1e-10, maxit=10000))
    while(any(resu$par!=resu0$par)){
        resu0 <- resu
        resu <- optim(par=resu0$par, fn=fn, control=list(factr = 1e-10, maxit=10000))
    }
    resu}
myoptimbounds <- function(par, fn, lower, upper, maxit=100){
    optim(par, fn=fn,
          gr = function(x) pracma::grad(fn, x), 
          method = "L-BFGS-B",
          lower = lower, upper = upper,
          control = list(factr = 1e-10, maxit = maxit))
}
mutualinfo <- function(jointFreqs, equalstim=FALSE, base=2L){##in bits by default
    if(equalstim){
        jointFreqs <- (jointFreqs/rowSums(jointFreqs))/nrow(jointFreqs)
    }else{
        jointFreqs <- jointFreqs/sum(jointFreqs)
    }
    jointFreqs[is.na(jointFreqs)] <- 0
    jointFreqs <- jointFreqs *
        log2(jointFreqs/outer(rowSums(jointFreqs), colSums(jointFreqs)))
    sum(jointFreqs[is.finite(jointFreqs)])/log2(base)
}
mutualinfoX <- function(jointFreqs, pstim1, base=2L){##in bits by default
        jointFreqs <- (jointFreqs/rowSums(jointFreqs))*c(1-pstim1, pstim1) # conditionals p(spike|stim) * p(stim)
        jointFreqs[is.na(jointFreqs)] <- 0
    jointFreqs <- jointFreqs *
        log2(jointFreqs/outer(rowSums(jointFreqs), colSums(jointFreqs)))
    sum(jointFreqs[is.finite(jointFreqs)])/log2(base)
}
capacity <- function(jointFreqs, base=2L){
    fn <- function(x){-mutualinfoX(jointFreqs, x)}
    -myoptimbounds(0.5, fn, 0, 1)$value/log2(base)
}
entropy <- function(freqs, base=2L){
    freqs <- cbind(freqs)
    freqs <- t(t(freqs)/colSums(freqs, na.rm=T))
    colSums(freqs*log2(1/freqs), na.rm=T)/log2(base)
}
## function to normalize absolute frequencies
normalize <- function(freqs){freqs/sum(freqs, na.rm=T)}
normalizerows <- function(freqs){freqs/rowSums(freqs, na.rm=T)}
normalizecols <- function(freqs){t(t(freqs)/colSums(freqs, na.rm=T))}

equalstim <- FALSE
set.seed(147)
binwidthms <- 500
longrunDataFile  <- 'BinarizedStimulus_SpikeCounts_dt=500ms.csv'
#sampleIndexFile  <- 'index_mat_80.csv'
T <- 16 ## prior weight
priorMeanSpikes <- 5*binwidthms/1000 # 5Hz
##
## load full recording
longrunData  <- as.data.table(t(read.csv(longrunDataFile,header=FALSE,sep=',')))
colnames(longrunData) <- c('stimulus','nspikes')
## longrunData$stimulus_lag1 <- NA
## longrunData$stimulus_lag1[2:nrow(longrunData)] <- 1L*longrunData$stimulus[2:nrow(longrunData)] + 2L*longrunData$stimulus[1:(nrow(longrunData)-1)]
##
stimuli <- sort(unique(longrunData[2:nrow(longrunData),stimulus],na.rm=T))
nStimuli <- length(stimuli)
######### shuffle
## longrunData <- longrunData[c(1,sample(2:nrow(longrunData)))] #shuffle
#########
## frequencies of full recording
maxSpikes <- 99#round(max(longrunData$nspikes,na.rm=T)*1.1)
maxSpikes1 <- maxSpikes + 1
#priorBaseDistr <- normalize(dpois(x=0:maxSpikes, lambda=priorMeanSpikes, log=FALSE))
priorBaseDistr <- normalize(dgeom(x=0:maxSpikes, prob=1/(priorMeanSpikes+1), log=FALSE))
## priorBaseDistr <- normalize(rep(1,maxSpikes1))
longrunFreqs <- t(sapply(stimuli, function(stim){
    tabulate(longrunData[stimulus==stim,nspikes]+1L, nbins=maxSpikes1)
}))
## Autocorrelation
pdff(paste0('autocorr_bin',binwidthms,'ms'))
lags <- 0:50
tplot(x=lags,y=cbind(coda::autocorr(coda::as.mcmc(longrunData$nspikes),lags),coda::autocorr(coda::as.mcmc(longrunData$stimulus),lags)),ylim=c(min(0,cbind(coda::autocorr(coda::as.mcmc(longrunData$nspikes),lags),coda::autocorr(coda::as.mcmc(longrunData$stimulus),lags))),NA),xlab='bins',ylab='autocorrelation',lwd=4)
legend('topright',legend=c('n. spikes','stimulus'),col=1:2,lty=1:2,bty='n',cex=2,lwd=3)
dev.off()
##
##priorBaseDistr <- normalize(c(colSums(longrunFreqs))+0.1)
rownames(longrunFreqs) <- nameStimulus <- paste0('stimulus',stimuli)
colnames(longrunFreqs) <- nameNspikes <- paste0('nspikes',0:maxSpikes)
longrunMI <- c(bit=mutualinfo(longrunFreqs, equalstim=equalstim))
longrunC <- capacity(longrunFreqs)

## addspikes <- longrunData$nspikes[seq(1,nrow(longrunData)-1,by=2)] +
##     longrunData$nspikes[seq(1,nrow(longrunData)-1,by=2)+1]
## addstims <- longrunData$stimulus[seq(1,nrow(longrunData)-1,by=2)]*2L +
##     longrunData$stimulus[seq(1,nrow(longrunData)-1,by=2)+1]
## longrunData2 <- data.table(nspikes=addspikes, stimulus=addstims)
## Autocorrelation
## pdff('autocorr_2bintogether')
## lags <- 0:50
## tplot(x=lags,y=cbind(coda::autocorr(coda::as.mcmc(longrunData2$nspikes),lags),coda::autocorr(coda::as.mcmc(longrunData2$stimulus),lags)),ylim=c(0,NA),xlab='bin (80 ms)',ylab='autocorrelation',lwd=4)
## legend('topright',legend=c('n. spikes','stimulus'),col=1:2,lty=1:2,bty='n',cex=2,lwd=3)
## dev.off()


nDraws <- 2^12
nPlotSamples <- 32
maxX <- maxSpikes
priorAlphas <- NULL
for(i in 1:nStimuli){priorAlphas <- rbind(priorAlphas, priorBaseDistr)}
dimnames(priorAlphas) <- list(nameStimulus, nameNspikes)
T <- 32 ## prior weight
startr <- 2 ## 2 or 0 for sequence with max uniformity of stimuli
startbin <- 1L
##
##pdff(paste0('MIhistogram_4stimEQ_prior',T,'_start',(if(startr==2){2}else{'OPT'})))
set.seed(149)
pdff(paste0('newMIhisto_bin',binwidthms,'ms_prior',T,'_startbin',startbin))
for(lsample in c(40*(1:(floor(nrow(longrunData)/40)-1)), nrow(longrunData))){
    ## if(FALSE){
    ##     entseq <- foreach(i=1L:(nrow(longrunData)-lsample+2L), .combine=c, .inorder=TRUE)%dopar%{
    ##         entropy(
    ##             tabulate(longrunData[i-1L+(1L:lsample),stimulus]+1L, nbins=nStimuli)
    ##         )
    ##     }
    ##     entorder <- order(entseq, decreasing=T)+1L
    ##     startbin <- entorder[1]
    ## }else{
    ##     startbin <- 1
    ## }
    ##
    sampleData <- longrunData[startbin-1+(1:lsample),]
    sampleFreqs <- t(sapply(stimuli, function(stim){
        tabulate(sampleData[stimulus==stim,nspikes]+1L, nbins=maxSpikes1)
    }))
    dimnames(sampleFreqs) <- dimnames(longrunFreqs)
    sampleMI <- c(bit=mutualinfo(sampleFreqs))
    sampleC <- c(bit=capacity(sampleFreqs))
    ##
    ## Alphas for Dirichlet distribution
    dAlphas <- sampleFreqs + T*priorAlphas/2
    ## Generate samples
    ## if(TRUE){
    ##     mcmcrun <- NULL
    ##     for(i in 1:nStimuli){mcmcrun <- cbind(mcmcrun, t(LaplacesDemon::rdirichlet(n=nDraws, alpha=dAlphas[i,])))}
    ##     dim(mcmcrun) <- c(maxSpikes1, nDraws, nStimuli)
    ##     mcmcrun <- aperm(mcmcrun, c(2,3,1))
    ## }else{
        mcmcrun <- LaplacesDemon::rdirichlet(n=nDraws, alpha=c(dAlphas))
        dim(mcmcrun) <- c(nDraws, nStimuli, maxSpikes1)
##        }
    dimnames(mcmcrun) <- list(NULL, nameStimulus, nameNspikes)
    postMISamples <- apply(mcmcrun, 1, function(x){mutualinfo(x)})
##    postCSamples <- apply(mcmcrun, 1, function(x){capacity(x)})
    postCSamples <- foreach(i=1:nDraws,.combine=c)%dopar%{capacity(mcmcrun[i,,])}
    ##
    postMIDistr <- thist(postMISamples)
    postCDistr <- thist(postCSamples)
    ##
    tplot(x=list(postMIDistr$mids, postCDistr$mids),y=list(postMIDistr$density, postCDistr$density),
          xlim=c(0,max(postMIDistr$breaks,postCDistr$breaks,longrunMI,longrunC,sampleMI,sampleC,0.75)),
          xlab='Sh', ylab='probability density',
          main=paste0(binwidthms,' ms; ',lsample, ' observed bins; start observation at bin ', startbin, '; prior weight: ', T),
          cex.main=1.25, col=c(1,6), lwd=2)
    tplot(x=list(longrunMI, longrunC, sampleMI, sampleC),
          y=list(-0.01*max(c(postMIDistr$density, postCDistr$density))),
    type='p', cex=2, pch=c(16,1,17,2), col=c(1,6,5,2), lty=c(1,2,1,2), add=T)
    ##    
    ## abline(v=longrunMI, col=3, lty=1, lwd=3)
    ## abline(v=longrunC, col=3, lty=2, lwd=3)
    ## abline(v=sampleMI, col=4, lty=1, lwd=2)
    ## abline(v=sampleC, col=4, lty=2, lwd=2)
    postMIQuantiles <- quant(x=postMISamples, probs=c(1,8,15)/16)
    postCQuantiles <- quant(x=postCSamples, probs=c(1,8,15)/16)
    ## abline(v=postMIQuantiles, col=c(5,5,5), lty=c(4,4,4), lwd=c(3,3,3))
    ## abline(v=postCQuantiles, col=c(2,2,2), lty=c(4,4,4), lwd=c(3,3,3))
    legend('topleft', legend=c(
                           'forecast long-run MI',
                           'forecast long-run capacity',
                           'true long-run MI',
                           'true long-run capacity',
                           'sample MI',
                           'sample capacity'),
           col=c(1,6,1,6,5,2),
           lty=c(1,2,NA,NA,NA,NA),
           pch=c(NA, NA, 16,1, 17,2),
           lwd=3, bty='n', cex=1)
    ##
    subsam <- sample(1:nrow(mcmcrun),size=nPlotSamples)
    ylim=c(-1,1)*max(mcmcrun[subsam,,])
    tplot(x=0:maxX, y=t(normalizerows(mcmcrun[subsam,1,])),type='l',lwd=0.5,lty=1,alpha=0.5,col=1,ylim=ylim)
    tplot(x=0:maxX, y=-t(normalizerows(mcmcrun[subsam,2,])),type='l',lwd=0.5,lty=1,alpha=0.5,col=1, add=T)
    tplot(x=0:maxX, y=t(c(1,-1)*normalizerows(longrunFreqs)), type='l', lwd=3, lty=1, col=c(3), add=T)
    tplot(x=0:maxX, y=t(c(1,-1)*normalizerows(sampleFreqs)), type='l', lwd=4, lty=2, col=c(4), add=T)
#    tplot(x=0:maxX, y=t(c(1,-1)*normalizerows(priorAlphas)), type='l', lwd=3, lty=3, col=7, add=T)
}
dev.off()



###################
#### OLD STUFF ####
###################

condfreqSamples <- t(apply(mcmcrun,1,normalizerows))
dim(condfreqSamples) <- c(nDraws, nStimuli, maxSpikes1)
dimnames(condfreqSamples) <- dimnames(mcmcrun)
    ##

    pdff(paste0('DDplots_samples',nSamples,'_chunk',chunk))
    ##
    matplot(x=0:maxSpikes, y=t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),1,]),
            type='l', lty=1, lwd=2, col=paste0(mypurpleblue,'22'), ylim=c(-1,1),  xlim=c(0,maxX),
            xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
    for(i in 2:nStimuli){
        matplot(x=0:maxSpikes, y=-t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),i,]),
                type='l', lty=1, lwd=1, col=paste0(myredpurple,'22'),
                add=TRUE)
    }
    matplot(x=0:maxSpikes, y=normalize(longrunFreqs[1,]),
            type='l', lty=1, lwd=2, col='black',
            add=TRUE)
    matplot(x=0:maxSpikes, y=-normalize(longrunFreqs[2,]),
            type='l', lty=1, lwd=2, col='black',
            add=TRUE)
    if(pflag==0){matplot(x=0:maxSpikes, y=normalize(sampleFreqs[1,]),
                         type='l', lty=2, lwd=5, col=myyellow,
                         add=TRUE)
                         matplot(x=0:maxSpikes, y=-normalize(sampleFreqs[2,]),
                                 type='l', lty=2, lwd=5, col=myyellow,
                                 add=TRUE)}
    title(paste0('(',nSamples,' data samples,',
                 ' chunk ', chunk,
                 ', prior weight = ', T, ')',
                 '\n superdistr ',chunk), cex.main=2)
    legend('topright',c('long-run freqs','sample freqs'),lty=c(1,2),lwd=c(2,5),col=c('black',myyellow),cex=1.5)
    ## Posterior MI
    matplot(x=postMIDistr$mids, y=postMIDistr$density,
            type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
            xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
    for(q in postMIQuantiles){
        matlines(x=rep(q,2),y=c(-1,1/2)*max(postMIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=4, lwd=6, col=myyellow)
    matlines(x=rep(longrunMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myredpurple)
    title('posterior MI distr', cex.main=2)
    legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myyellow,myredpurple),lwd=4,cex=1.5)
    dev.off()
    ##
    NULL
}


## frequencies of sample
gc()
nores <- foreach(chunk=0:2, .inorder=F, .packages=c('data.table'))%dorng%{
    pflag <- 0
    if(chunk==0){chunk <- 1
        pflag <- 1}
    chunkIndices <- as.matrix(read.csv(sampleIndexFile,header=FALSE,sep=','))[chunk,]
    sampleData <- longrunData[chunkIndices,]
    ##print(str(sampleData))
    sampleFreqs <- foreach(stim=stimuli, .combine=rbind)%do%{
        tabulate(sampleData[stimulus==stim,nspikes]+1, nbins=maxSpikes1)
    } 
    dimnames(sampleFreqs) <- dimnames(longrunFreqs)
    nSamples <- sum(sampleFreqs)
    if(pflag==0){
        sampleMI <- c(bit=mutualinfo(sampleFreqs))
    } else {
        sampleMI <- c(bit=-2)
        chunk <- 0}
    ##
    ##
    ## Alphas for Dirichlet distribution
    ##
    priorAlphas <- normalize(matrix(rep(priorBaseDistr, each=nstimuli),nrow=nstimuli))
    dimnames(priorAlphas) <- list(nameStimulus, nameNspikes)
    dAlphas <- T*priorAlphas + if(pflag==0){sampleFreqs}else{0}
    ## Generate samples
    set.seed(147+chunk)
    mcmcrun <- rdirichlet(n=nDraws, alpha=c(dAlphas))
    dim(mcmcrun) <- c(nDraws,nstimuli,maxSpikes1)
    dimnames(mcmcrun) <- list(NULL, nameStimulus, nameNspikes)
    postMISamples <- apply(mcmcrun,1,mutualinfo)
    postMIDistr <- hist(postMISamples, breaks=seq(0,1,by=0.02), plot=F)
    postMIQuantiles <- quantile(x=postMISamples, probs=c(0.025,0.5,0.975))
    condfreqSamples <- t(apply(mcmcrun,1,normalizerows))
    dim(condfreqSamples) <- c(nDraws, nStimuli, maxSpikes1)
    dimnames(condfreqSamples) <- dimnames(mcmcrun)
    ##
    ##
    ## Plot draws
    pdff(paste0('DDplots_samples',nSamples,'_chunk',chunk))
    ##
    matplot(x=0:maxSpikes, y=t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),1,]),
            type='l', lty=1, lwd=2, col=paste0(mypurpleblue,'22'), ylim=c(-1,1),  xlim=c(0,maxX),
            xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
    for(i in 2:nStimuli){
        matplot(x=0:maxSpikes, y=-t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),i,]),
                type='l', lty=1, lwd=1, col=paste0(myredpurple,'22'),
                add=TRUE)
    }
    matplot(x=0:maxSpikes, y=normalize(longrunFreqs[1,]),
            type='l', lty=1, lwd=2, col='black',
            add=TRUE)
    matplot(x=0:maxSpikes, y=-normalize(longrunFreqs[2,]),
            type='l', lty=1, lwd=2, col='black',
            add=TRUE)
    if(pflag==0){matplot(x=0:maxSpikes, y=normalize(sampleFreqs[1,]),
                         type='l', lty=2, lwd=5, col=myyellow,
                         add=TRUE)
                         matplot(x=0:maxSpikes, y=-normalize(sampleFreqs[2,]),
                                 type='l', lty=2, lwd=5, col=myyellow,
                                 add=TRUE)}
    title(paste0('(',nSamples,' data samples,',
                 ' chunk ', chunk,
                 ', prior weight = ', T, ')',
                 '\n superdistr ',chunk), cex.main=2)
    legend('topright',c('long-run freqs','sample freqs'),lty=c(1,2),lwd=c(2,5),col=c('black',myyellow),cex=1.5)
    ## Posterior MI
    matplot(x=postMIDistr$mids, y=postMIDistr$density,
            type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
            xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
    for(q in postMIQuantiles){
        matlines(x=rep(q,2),y=c(-1,1/2)*max(postMIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=4, lwd=6, col=myyellow)
    matlines(x=rep(longrunMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myredpurple)
    title('posterior MI distr', cex.main=2)
    legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myyellow,myredpurple),lwd=4,cex=1.5)
    dev.off()
    ##
    NULL
}
