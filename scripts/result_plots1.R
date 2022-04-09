## Author: PGL  Porta Mana
## Created: 2022-01-12T14:51:16+0100
## Last-Updated: 2022-04-10T00:19:28+0200
################
## Relation between softmax & probability
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
## library('LaplacesDemon')
## library('extraDistr')
## library('mvtnorm')
## options(bitmapType='cairo')
## pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
## pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
library('nimble')
#### End custom setup ####

## Bernoulli distribution
dbernoulli <- function(x, prob, log=FALSE){
    if(log){
        out <- x*log(prob) + (1-x)*log(1-prob)
        out[is.na(out)] <- 0
    }else{
        out <- x*prob + (1-x)*(1-prob)
    }
    out
}

maincov <- 'class'
source('functions_mcmc.R')
dirname <- '_testmcmc1_-V2-D4096-K16-I4096'
frequenciesfile <- '_mcsamples-Rtestmcmc1_2-V2-D4096-K16-I4096.rds'
##
datafile <- 'softmaxdata_test.csv'
varinfofile <- 'variate_info.csv'
##
variateinfo <- fread(varinfofile, sep=',')
covNames <- variateinfo$variate
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- covNames
realCovs <- covNames[covTypes=='double']
integerCovs <- covNames[covTypes=='integer']
binaryCovs <- covNames[covTypes=='binary']
covNames <- c(realCovs, integerCovs, binaryCovs)
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)
##
frequenciesfile <- paste0(dirname,'/',frequenciesfile)
outfile <- paste0(dirname,'/','results.txt')
parmList <- mcsamples2parmlist(readRDS(frequenciesfile))
nclusters <- ncol(parmList$q)
nFsamples <- nrow(parmList$q)
##
otherCovs <- setdiff(covNames, maincov)
##
alldata <- fread(datafile, sep=',')
alldata <- alldata[,..covNames]
#alldata <- alldata[Usage_ == 'train']

## probc0 <- samplesF(Y=matrix(0,nrow=1,dimnames=list(NULL,maincov)), X=data.matrix(alldata[2,..otherCovs]), parmList=parmList, inorder=T)
## ##
## probc1 <- samplesF(Y=matrix(1,nrow=1,dimnames=list(NULL,maincov)), X=data.matrix(alldata[2,..otherCovs]), parmList=parmList, inorder=T)
## ##
## probj0 <- samplesF(X=NULL, Y=cbind(matrix(0,nrow=1,dimnames=list(NULL,maincov)),data.matrix(alldata[2,..otherCovs])), parmList=parmList, inorder=T)
## probj1 <- samplesF(X=NULL, Y=cbind(matrix(1,nrow=1,dimnames=list(NULL,maincov)),data.matrix(alldata[2,..otherCovs])), parmList=parmList, inorder=T)
## ##
## probx <- samplesF(X=NULL, Y=data.matrix(alldata[2,..otherCovs]), parmList=parmList, inorder=T)

smgrid <- matrix(seq(1e-6, 1-1e-6, length.out=256), ncol=1, dimnames=list(NULL,'softmax'))
lsmgrid <- qlogis(smgrid)
colnames(lsmgrid) <- 'logitsoftmax'
##
probf <- samplesF(Y=cbind(class=0), X=lsmgrid, parmList=parmList, inorder=F)
##
qprobf <- apply(probf,1,function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
##
pdff(paste0(dirname,'/','prob_vs_softmax2'))
ylim <- c(0,1)
xlim <- c(0,1)
##
tplot(x=smgrid, y=qprobf[1,], yticks=NULL, xlim=xlim,
      lwd=3, alpha=0.25, ylim=ylim, cex.axis=1.25, cex.lab=1.5,
      xlab='softmax output 0', ylab='probability of class 0')
polygon(x=c(smgrid,rev(smgrid)), y=c(qprobf[2,], rev(qprobf[3,])), col=paste0(palette()[1],'40'), border=NA)
    ## legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients who will have ',diseasenames), '87.5% uncertainty'),
    ##        col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
    ##        )
    ## for(i in 1:2){
    ##     histo <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=32)
    ##     tplot(x=histo$breaks,y=histo$density,col=i,border=NA,alpha=0.75,xgrid=F,ygrid=F,add=T)
    ## }
dev.off()

ndata <- 4096L
schoices <- 1L*(!(alldata$logitsoftmax > 0))
probd <- sapply(0:1, function(xclass){sapply(0:1, function(xsch){
    sum(alldata$class[1:ndata]==xclass & schoices[1:ndata]==xsch)
})})
dimnames(probd) <- list(paste0('schoice.',0:1),paste0('true.',0:1))
##
cprob <- probd/rowSums(probd)
cprob0 <- cprob[schoices+1,1]
names(cprob0) <- NULL

confmatrix <- sapply(0:1, function(xclass){sapply(0:1, function(xsch){
    sum(alldata$class==xclass & schoices==xsch)
})})
dimnames(confmatrix) <- list(paste0('schoice.',0:1), paste0('true.',0:1))



testdata <- tail(alldata,n=8192)
transfsm <- samplesF(Y=cbind(class=0), X=cbind(logitsoftmax=testdata$logitsoftmax), parmList=parmList, inorder=T)
prob0 <- rowMeans(transfsm)
testdata <- cbind(testdata, data.table(prob=prob0))

decidevaluate <- function(truevalues, probs0, umatrix, normalize=F, shift=F, average=F){
    dimnames(umatrix) <- list(paste0('choice.',0:1), paste0('true.',0:1))
    if(shift){ umatrix <- umatrix - min(umatrix) }
    if(normalize){ umatrix <- umatrix/max(abs(umatrix)) }
    if(normalize | shift){
        print('rescaled utility matrix:')
        print(umatrix)
    }
    print(paste0('probability threshold: ',
(umatrix[2,2]-umatrix[1,2])/(umatrix[1,1]+umatrix[2,2]-umatrix[1,2]-umatrix[2,1])
                 ))
    ##
    exputilities <- umatrix %*% rbind(probs0, 1-probs0)
    decmatrix <- rbind(c(
        sum(truevalues==0 & exputilities[1,]>exputilities[2,]) +
        sum(truevalues==0 & exputilities[1,]==exputilities[2,])/2 ,
        sum(truevalues==1 & exputilities[1,]>exputilities[2,]) +
        sum(truevalues==1 & exputilities[1,]==exputilities[2,])/2
        ), c(
        sum(truevalues==0 & exputilities[1,]<exputilities[2,]) +
        sum(truevalues==0 & exputilities[1,]==exputilities[2,])/2 ,
        sum(truevalues==1 & exputilities[1,]<exputilities[2,]) +
        sum(truevalues==1 & exputilities[1,]==exputilities[2,])/2
        ))
    dimnames(decmatrix) <- dimnames(umatrix)
    print('decision matrix:')
    print(decmatrix)
    ##
    res <- sum(decmatrix * umatrix)
    if(average){res <- res/length(truevalues)}
    print(paste0((if(average){'average'}else{'total'}), ' utility:'))
    print(res)
    res
}

decidevaluate(testdata$class, plogis(testdata$logitsoftmax), diag(2),average=F)
decidevaluate(testdata$class, testdata$prob, diag(2),average=F)

um <- rbind(c(1,-99),c(0,1))
decidevaluate(testdata$class, plogis(testdata$logitsoftmax), um,average=F)
decidevaluate(testdata$class, testdata$prob, um,average=F)

um <- rbind(c(1,0),c(-99,1))
decidevaluate(testdata$class, plogis(testdata$logitsoftmax), um,average=F)
decidevaluate(testdata$class, testdata$prob, um,average=F)

um <- rbind(c(0.01,0),c(0,1))
decidevaluate(testdata$class, plogis(testdata$logitsoftmax), um,average=F)
decidevaluate(testdata$class, testdata$prob, um,average=F)

um <- rbind(c(1,0),c(0,0.01))
decidevaluate(testdata$class, plogis(testdata$logitsoftmax), um,average=F)
decidevaluate(testdata$class, testdata$prob, um,average=F)


for(um in list(diag(2),
               rbind(c(1,-999),c(0,1)),
               rbind(c(1,0),c(-999,1)),
               rbind(c(0.001,0),c(0,1)),
               rbind(c(1,0),c(0,0.001))
               )){
    print('utility matrix:')
    print(um)
    resu <- decidevaluate(alldata$class, cprob0, um, average=F)
    print('softmax utility:')
    print(c(sum(um * confmatrix)))
}


                                        #um <- diag(2)

um <- 

decidevaluate(testdata$class, plogis(testdata$logitsoftmax), um,average=F)
decidevaluate(testdata$class, testdata$prob, um,average=F)

um <- 
decidevaluate(testdata$class, plogis(testdata$logitsoftmax), um,average=F)
decidevaluate(testdata$class, testdata$prob, um,average=F)

um <- 
decidevaluate(testdata$class, plogis(testdata$logitsoftmax), um,average=F)
decidevaluate(testdata$class, testdata$prob, um,average=F)

um <- 
decidevaluate(testdata$class, plogis(testdata$logitsoftmax), um,average=F)
decidevaluate(testdata$class, testdata$prob, um,average=F)




## formula for threshold
## u00 Q + u01 (1-Q) > u10 Q + u11 (1-Q)
## Q (u00-u01-u10+u11) > u11-u01
## Q > (u11-u01)/(u00+u11-u10-u01)







#################
## Old code below


tplot(x=smgrid, y=rowMeans(probf))

grids <- foreach(acov=covNames)%do%{
    rg <- range(alldata[[acov]], na.rm=T)
    if(acov %in% realCovs){
        rg <- rg+c(-1,1)*IQR(alldata[[acov]],type=8,na.rm=T)/2
        Xgrid <- seq(rg[1], rg[2], length.out=256)
    }else if(acov %in% integerCovs){
            rg <- round(c((covMins[acov]+3*rg[1])/4, (covMaxs[acov]+3*rg[2])/4))
            Xgrid <- rg[1]:rg[2]
    }else{Xgrid <- 0:1}
    matrix(Xgrid, ncol=1, dimnames=list(NULL,acov))
}
names(grids) <- covNames
##
xcond <- matrix(0:1,ncol=2,dimnames=list(NULL,rep(maincov,2)))

## Frequencies of each feature given AD state
distsFA <- foreach(acov=otherCovs)%do%{
    dists <- rbind(samplesF(Y=grids[[acov]], X=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(Y=grids[[acov]], X=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, diseasenames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsFA) <- otherCovs

## quantiles
qdistsFA <- foreach(acov=otherCovs)%do%{
    apply(distsFA[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qdistsFA) <- otherCovs

#allMI <- readRDS(paste0(dirname,'/allMI.rds'))
orderc <- order(apply(allMI[otherCovs,],1,median),decreasing=T)
svnames <- sapply(covNames,function(acov){gsub('([^_]+)_.*', '\\1', acov)})

## plot of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','allplots_features_given_AD'))
par(mfrow=c(3,4), mai=c(0,0,0,0),oma=c(0,0,0,0))
#par(mfrow=c(3,4), mai=c(0,0,0,0),oma=c(0,0,0,0),mar=c(0,0,0,0))
for(acov in otherCovs[orderc]){
    agrid <- grids[[acov]]
    ymax <- quant(apply(qdistsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
    ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    for(i in 1:2){
        tplot(x=agrid, y=qdistsFA[[acov]][1,,i], yticks=NULL, xlim=xlim,
              col=i, lty=i, lwd=3, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,ylabels=F,cex.axis=1.25,cex.lab=1.5,ly=0,
              xlab=svnames[acov], ylab=NA, add=(i==2))
        polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
    }
    ## legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients who will have ',diseasenames), '87.5% uncertainty'),
    ##        col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
    ##        )
    ## for(i in 1:2){
    ##     histo <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=32)
    ##     tplot(x=histo$breaks,y=histo$density,col=i,border=NA,alpha=0.75,xgrid=F,ygrid=F,add=T)
    ## }
}
dev.off()
## plot of samples of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','exampleplotsamples_features_given_AD'))
for(acov in otherCovs[orderc[1]]){
    agrid <- grids[[acov]]
    ymax <- quant(apply(distsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
    ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    subsam <- seq(1,dim(distsFA[[acov]])[1], length.out=128)
    tplot(x=agrid, y=matrix(
                       rbind(t(distsFA[[acov]][subsam,,1]),t(distsFA[[acov]][subsam,,2])),
                       nrow=length(agrid)),
          yticks=NULL, xlim=xlim,
              col=c(5,2), lty=1, lwd=1, alpha=0.75, ylim=ylim, xticks=xticks, xlabels=xlabels,
              xlab=svnames[acov], ylab='population frequency')
}
dev.off()
## plot of samples of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','exampleplot_features_given_AD'))
for(acov in otherCovs[orderc[1]]){
    agrid <- grids[[acov]]
    ##ymax <- quant(apply(qdistsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
    ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    for(i in 1:2){
        tplot(x=agrid, y=qdistsFA[[acov]][1,,i], yticks=NULL, xlim=xlim,
              col=i, lty=i, lwd=4, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,
              xlab=svnames[acov], ylab='population frequency', add=(i==2))
        polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
    }
    ## legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients who will have ',diseasenames), '87.5% uncertainty'),
    ##        col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
    ##        )
    ## for(i in 1:2){
    ##     histo <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=32)
    ##     tplot(x=histo$breaks,y=histo$density,col=i,border=NA,alpha=0.75,xgrid=F,ygrid=F,add=T)
    ## }
}
dev.off()



## probability of AD state given features, via Bayes's theorem
bayesAF <- foreach(acov=otherCovs)%do%{
    dist <- distsFA[[acov]]
    zz <- dist[,,'AD']+dist[,,'MCI']
    dist[,,'AD'] <- dist[,,'AD']/zz
    dist[,,'MCI'] <- dist[,,'MCI']/zz
    dist
}
names(bayesAF) <- otherCovs

## quantiles
qbayesAF <- foreach(acov=otherCovs)%do%{
    apply(bayesAF[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qbayesAF) <- otherCovs

## frequencies of AD state given features
distsAF <- foreach(acov=setdiff(covNames, maincov))%do%{
    dists <- rbind(samplesF(X=grids[[acov]], Y=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(X=grids[[acov]], Y=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, diseasenames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsAF) <- otherCovs


## quantiles
qdistsAF <- foreach(acov=otherCovs)%do%{
    apply(distsAF[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qdistsAF) <- otherCovs

## plot of frequencies of AD state given features, f(AD|F)
pdff(paste0(dirname,'/','allplots_predictAD'))
par(mfrow=c(3,4), mai=c(0,0,0,0),oma=c(0,0,0,0))
for(iacov in 1:length(otherCovs[orderc])){acov <- otherCovs[orderc][iacov]
    agrid <- grids[[acov]]
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
ylim <- c(0,1)
    xlim <- range(agrid)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    tplot(x=c(-100,100), y=c(-100,100),type='p',
          col=2, lty=1, lwd=3, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          yticks=(0:4)/4, ylabels=(if((iacov-1)%%4==0){c('0%','25%','50%','75%','100%')}else{F}),cex.axis=1.25,cex.lab=1.5,ly=(if((iacov-1)%%4==0){2.5}else{0}),
          xlab=svnames[acov], ylab=NA)#(if((iacov-1)%%4==0){'probability of AD onset'}else{NA}))
    polygon(x=c(xlim,rev(xlim)),y=c(0.5,0.5,1,1),col=paste0(palette()[2],'22'),border=NA)
    polygon(x=c(xlim,rev(xlim)),y=c(0,0,0.5,0.5),col=paste0(palette()[1],'22'),border=NA)
    abline(h=0.5, lty=2, lwd=1, col=4)
    tplot(x=agrid, y=qdistsAF[[acov]][1,,'AD'],
          col=darkgrey, lty=1, lwd=3, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          yticks=(0:4)/4, ylabels=(if((iacov-1)%%4==0){c('0%','25%','50%','75%','100%')}else{F}),cex.axis=1.25,cex.lab=1.5,ly=(if((iacov-1)%%4==0){2.5}else{0}),
          xlab=svnames[acov], ylab=NA,add=T)#(if((iacov-1)%%4==0){'probability of AD onset'}else{NA}))
    #mtext('probability of AD onset', side=2, outer=T, at=0.5,cex=1)
    polygon(x=c(agrid,rev(agrid)), y=c(qdistsAF[[acov]][2,,'AD'], rev(qdistsAF[[acov]][3,,'AD'])), col=paste0(darkgrey,'80'), border=NA)
## legend(x=agrid[1], y=ylim[2]*1, legend=c('87.5% uncertainty on the probability'),
##        col=palette()[c(2)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
##                        )
}
dev.off()



###################################
#### Graphs


## plot of frequencies of features given AD state f(F|AD)
withprob <- !TRUE
pdff(paste0(dirname,'/','histograms_and_probs-withprobs',withprob))
par(mfrow=c(3,4), mai=c(0,0,0,0),oma=c(0,0,0,0))
for(acov in otherCovs[orderc]){
    agrid <- grids[[acov]]
    histo <- list()
##    breaks <- if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=24)}else{seq(min(agrid)-0.5,max(agrid)+0.5,by=1)}
    breaks <- if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=32)}else{seq(min(agrid)-0.5,max(agrid)+0.5,length.out=32)}
    ## breaks <- (if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=32)}else{NULL})
    ymax <- -Inf
    for(i in 1:2){
        histo[[i]] <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])), n=breaks)
        ymax <- max(ymax, histo[[i]]$counts)
            }
    ylim <- c(-ymax*1.1,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    freqs <- histo[[2]]$counts/(histo[[1]]$counts+histo[[2]]$counts)
    maxf <- 1#max(freqs,na.rm=T)
    freqs <- freqs/maxf*ymax
    yticks <- c((((-4):0)/4/maxf*ymax)-0.1/maxf*ymax,pretty(c(0,ymax),5))
    ylabels <- c('100%','75%','50%','25%','0%',sprintf('%.7g',c(yticks[yticks>=0])))
    for(i in 1:2){
        ## histo <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=32)}else{'i'})
        tplot(x=histo[[i]]$breaks,y=histo[[i]]$counts,col=i,alpha=0.75,border=i,lwd=0.75, cex.lab=1.25,ly=2.5,lty=i,
              xlim=xlim,
              ylim=ylim, xticks=xticks, xlabels=xlabels, cex.axis=1,
              yticks=yticks, ylabels=ylabels,
              xlab=svnames[acov], ylab=NA, add=(i==2))
        ## polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
    }
    tplot(x=histo[[1]]$mids[!is.na(freqs)],y=-freqs[!is.na(freqs)]-0.1/maxf*ymax,col=4,xgrid=F,ygrid=F,add=T)
#    abline(h=-0.5/maxf*ymax,lty=3,lwd=3,col=4)
    ##
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
    ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
if(withprob){
    tplot(x=agrid, y=-qdistsAF[[acov]][1,,'AD']/maxf*ymax-0.1/maxf*ymax,
          col=7, lty=1, lwd=2, ylim=ylim, xlim=xlim, xticks=F,yticks=F,xgrid=F,ygrid=F,add=T)
    ## tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    ## if(!any(is.na(tpar))){
    ##     Ogrid <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
    ##     axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
    ## }
    polygon(x=c(agrid,rev(agrid)), y=-c(qdistsAF[[acov]][2,,'AD'], rev(qdistsAF[[acov]][3,,'AD']))/maxf*ymax-0.1/maxf*ymax, col=paste0(palette()[7],'80'), border=NA)
## legend('topleft', legend=c('87.5% uncertainty on the probability'),
##        col=palette()[c(7)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
    ##                        )
}
    ## tplot(x=c(data.matrix(alldata[[acov]])),y=rep(-0.02,nrow(alldata))/maxf*ymax,type='p',pch=15,cex=0.75,alpha=1-1/16,col='#000000',add=T,xgrid=F,ygrid=F,xticks=F,yticks=F)
    ## legend(x='topright', legend=c(paste0('',diseasenames)),
    ##        col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
    ##        )
}
dev.off()



## plot of frequencies of features given AD state f(F|AD)
withprob <- TRUE
pdff(paste0(dirname,'/','LRHHC-histograms_and_probs-withprobs',withprob))
for(acov in 'LRHHC_n_long_log'){
    agrid <- grids[[acov]]
    histo <- list()
##    breaks <- if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=24)}else{seq(min(agrid)-0.5,max(agrid)+0.5,by=1)}
    breaks <- if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=32)}else{seq(min(agrid)-0.5,max(agrid)+0.5,length.out=32)}
    ## breaks <- (if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=32)}else{NULL})
    ymax <- -Inf
    for(i in 1:2){
        histo[[i]] <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])), n=breaks)
        ymax <- max(ymax, histo[[i]]$counts)
            }
    ylim <- c(-ymax*1.1,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    freqs <- histo[[2]]$counts/(histo[[1]]$counts+histo[[2]]$counts)
    maxf <- 1#max(freqs,na.rm=T)
    freqs <- freqs/maxf*ymax
    yticks <- c((((-4):0)/4/maxf*ymax)-0.1/maxf*ymax,pretty(c(0,ymax),5))
    ylabels <- c('100%','75%','50%','25%','0%',sprintf('%.7g',c(yticks[yticks>=0])))
    for(i in 1:2){
        ## histo <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=32)}else{'i'})
        tplot(x=histo[[i]]$breaks,y=histo[[i]]$counts,col=i,alpha=0.75,border=i,lwd=0.75, cex.lab=1.25,ly=2.5,lty=i,
              xlim=xlim,
              ylim=ylim, xticks=xticks, xlabels=xlabels, cex.axis=1,
              yticks=yticks, ylabels=ylabels,
              xlab=svnames[acov], ylab=NA, add=(i==2))
        ## polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
    }
    tplot(x=histo[[1]]$mids[!is.na(freqs)],y=-freqs[!is.na(freqs)]-0.1/maxf*ymax,col=4,xgrid=F,ygrid=F,add=T)
#    abline(h=-0.5/maxf*ymax,lty=3,lwd=3,col=4)
    ##
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
    ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
if(withprob){
    tplot(x=agrid, y=-qdistsAF[[acov]][1,,'AD']/maxf*ymax-0.1/maxf*ymax,
          col=7, lty=1, lwd=2, ylim=ylim, xlim=xlim, xticks=F,yticks=F,xgrid=F,ygrid=F,add=T)
    ## tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    ## if(!any(is.na(tpar))){
    ##     Ogrid <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
    ##     axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
    ## }
    polygon(x=c(agrid,rev(agrid)), y=-c(qdistsAF[[acov]][2,,'AD'], rev(qdistsAF[[acov]][3,,'AD']))/maxf*ymax-0.1/maxf*ymax, col=paste0(palette()[7],'80'), border=NA)
## legend('topleft', legend=c('87.5% uncertainty on the probability'),
##        col=palette()[c(7)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
    ##                        )
}
    ## tplot(x=c(data.matrix(alldata[[acov]])),y=rep(-0.02,nrow(alldata))/maxf*ymax,type='p',pch=15,cex=0.75,alpha=1-1/16,col='#000000',add=T,xgrid=F,ygrid=F,xticks=F,yticks=F)
    ## legend(x='topright', legend=c(paste0('',diseasenames)),
    ##        col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
    ##        )
}
dev.off()




## plot of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','data_histogramsprobs_given_AD'),paper='a4')
for(acov in otherCovs){
    agrid <- grids[[acov]]
    histo <- list()
    breaks <- if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=24)}else{seq(min(agrid)-0.5,max(agrid)+0.5,by=1)}
    ymax <- -Inf
    for(i in 1:2){
        histo[[i]] <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=breaks)
        ymax <- max(ymax, histo[[i]]$counts)
            }
    ylim <- c(-ymax,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    freqs <- histo[[2]]$counts/(histo[[1]]$counts+histo[[2]]$counts)
    maxf <- 1#max(freqs,na.rm=T)
    freqs <- freqs/maxf*ymax
    yticks <- c(pretty(c(-1,0),10)/maxf*ymax,pretty(c(0,ymax),10))
    ylabels <- sprintf('%.7g',c(-yticks[yticks<0]*maxf/ymax, yticks[yticks>=0]))
    for(i in 1:2){
        ## histo <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=if(acov %in% realCovs){seq(min(agrid),max(agrid),length.out=32)}else{'i'})
        tplot(x=histo[[i]]$breaks,y=histo[[i]]$counts,col=i,alpha=0.75,
              xlim=xlim,
              ylim=ylim, xticks=xticks, xlabels=xlabels,
              yticks=yticks, ylabels=ylabels,
              xlab=acov, ylab='data counts', add=(i==2))
        ## polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
    }
    tplot(x=histo[[1]]$mids[!is.na(freqs)],y=-freqs[!is.na(freqs)],col=3,xgrid=F,ygrid=F,add=T)
    abline(h=-0.5/maxf*ymax,lty=3,lwd=3,col=4)
    ##
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
    ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    tplot(x=agrid, y=-qdistsAF[[acov]][1,,'AD']/maxf*ymax,
          col=7, lty=1, lwd=4, ylim=ylim, xlim=xlim, xticks=F,yticks=F,xgrid=F,ygrid=F,add=T)
    ## tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    ## if(!any(is.na(tpar))){
    ##     Ogrid <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
    ##     axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
    ## }
    polygon(x=c(agrid,rev(agrid)), y=-c(qdistsAF[[acov]][2,,'AD'], rev(qdistsAF[[acov]][3,,'AD']))/maxf*ymax, col=paste0(palette()[7],'80'), border=NA)
## legend('topleft', legend=c('87.5% uncertainty on the probability'),
##        col=palette()[c(7)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
    ##                        )
    tplot(x=c(data.matrix(alldata[[acov]])),y=rep(-0.02,nrow(alldata))/maxf*ymax,type='p',pch=15,cex=0.75,alpha=1-1/16,col='#000000',add=T,xgrid=F,ygrid=F,xticks=F,yticks=F)
    legend(x='topright', legend=c(paste0('',diseasenames)),
           col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
           )
}
dev.off()

##########################################################################


## plot of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','plots_features_given_AD'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ymax <- quant(apply(qdistsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
    ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    for(i in 1:2){
        tplot(x=agrid, y=qdistsFA[[acov]][1,,i], yticks=NULL, xlim=xlim,
              col=i, lty=i, lwd=4, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,
              xlab=acov, ylab='frequency of feature for patients with AD/MCI', add=(i==2))
        polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
    }
    legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients who will have ',diseasenames), '87.5% uncertainty'),
           col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
           )
    ## for(i in 1:2){
    ##     histo <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=32)
    ##     tplot(x=histo$breaks,y=histo$density,col=i,border=NA,alpha=0.75,xgrid=F,ygrid=F,add=T)
    ## }
}
dev.off()

## ## Data histogram plots for LRHHC
## tpar <- unlist(variateinfo[variate=='LRHHC_n_long_log',c('transfM','transfW')])
##         xlabels <- signif(pretty(exp(tpar['transfW']*hgrid + tpar['transfM']),n=10),2)
##         xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
## ## tplot(x=hgrid, y=lhist1$density/(lhist1$density+lhist0$density), ylim=c(0,NA),xticks=xticks,xlabels=xlabels)
## tplot(x=hgrid, y=lhist0$density, ylim=c(-1,NA),xticks=xticks,xlabels=xlabels,
##       xlab='LRHHC_n_long_log', ylab='P(AD)                                                                  frequencies')
## tplot(x=hgrid, y=lhist1$density, col=2, ylim=c(0,NA),xticks=xticks,xlabels=xlabels,add=T)
## tplot(x=hgrid, y=-lhist1$density/(lhist1$density+lhist0$density), col=3, ylim=c(0,NA),xticks=xticks,xlabels=xlabels,add=T)




## plot of samples of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','plotssamples_features_given_AD'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ymax <- quant(apply(distsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
    ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    subsam <- seq(1,dim(distsFA[[acov]])[1], length.out=128)
    tplot(x=agrid, y=matrix(
                       rbind(t(distsFA[[acov]][subsam,,1]),t(distsFA[[acov]][subsam,,2])),
                       nrow=length(agrid)),
          yticks=NULL, xlim=xlim,
              col=c(5,2), lty=1, lwd=1, alpha=0.75, ylim=ylim, xticks=xticks, xlabels=xlabels,
              xlab=acov, ylab='frequency of feature for patients with AD/MCI')
    ## for(i in 1:2){
    ##     tplot(x=agrid, y=qdistsFA[[acov]][1,,i], yticks=NULL, xlim=xlim,
    ##           col=c(1,6)[i], lty=i, lwd=3, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,
    ##           xlab=acov, ylab='frequency of feature for patients with AD/MCI', add=T)
    ## }
    legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distributions for patients who will have ',diseasenames)),
           col=palette()[c(1,2,7)], lty=c(1,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
           )
}
dev.off()

## ## plot of frequencies of features given AD state and gender f(F|AD&G)
## pdff('plots_features_given_ADG2')
## for(acov in names(distsFAG)){
##     agrid <- grids[[acov]]
##    ymax <- quant(apply(qdistsFAG[[acov]],2,function(x){quant(x,99/100)}),99/100)
## ##    ymax <- max(qdistsFAG[[acov]])
##     ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
##     xlim <- c(NA,NA)
##     tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
##     if(!any(is.na(tpar))){
##         xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
##         xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
##     }else{xticks <- NULL
##         xlabels <- TRUE}
##     if(acov %in% binaryCovs){
##         xticks <- 0:1
##         xlim <- c(-0.25,1.25)
##     }
##     tcols <- matrix(c(1,6,5,2),nrow=2)
##     for(i in 1:2){
##         for(j in 1:2){
##         tplot(x=agrid, y=qdistsFAG[[acov]][1,,i,j], yticks=NULL, xlim=xlim,
##               col=tcols[i,j], lty=i, lwd=4, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,
##               xlab=acov, ylab='frequency of feature given AD/MCI & gender', add=(i+j>2))
##         polygon(x=c(agrid,rev(agrid)), y=c(qdistsFAG[[acov]][2,,i,j], rev(qdistsFAG[[acov]][3,,i,j])), col=paste0(palette()[tcols[i,j]],'40'), border=NA)
##     }
##     }
##     legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients with ',diseasenames), '87.5% uncertainty'),
##            col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
##            )
##     legend(x=agrid[length(agrid)*4/5], y=ylim[2]*1.2, legend=c('darker: male','lighter: female'), bty='n', xpd=T, cex=1.25)
## }
## dev.off()


## plot of frequencies of AD state given features, using bayes f(F|AD)
pdff(paste0(dirname,'/','plots_predictAD_bayes'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
    ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    tplot(x=agrid, y=qbayesAF[[acov]][1,,'AD'],
          col=7, lty=1, lwd=4, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          xlab=acov, ylab='probability of AD')
    ## tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    ## if(!any(is.na(tpar))){
    ##     Ogrid <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
    ##     axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
    ## }
    polygon(x=c(agrid,rev(agrid)), y=c(qbayesAF[[acov]][2,,'AD'], rev(qbayesAF[[acov]][3,,'AD'])), col=paste0(palette()[7],'80'), border=NA)
    abline(h=0.5, lty=2, lwd=1, col=2)
legend('topleft', legend=c('87.5% uncertainty on the probability'),
       col=palette()[c(7)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
                       )
}
dev.off()

## plot of frequencies of AD state given features, f(AD|F)
pdff(paste0(dirname,'/','plots_predictAD'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    tplot(x=agrid, y=qdistsAF[[acov]][1,,'AD'],
          col=2, lty=1, lwd=4, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          xlab=acov, ylab='probability of AD')
    polygon(x=c(agrid,rev(agrid)), y=c(qdistsAF[[acov]][2,,'AD'], rev(qdistsAF[[acov]][3,,'AD'])), col=paste0(palette()[2],'80'), border=NA)
    abline(h=0.5, lty=2, lwd=1, col=2)
legend(x=agrid[1], y=ylim[2]*1, legend=c('87.5% uncertainty on the probability'),
       col=palette()[c(2)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
                       )
}
dev.off()

## plot of samples of frequencies of AD state given features, f(AD|F)
pdff(paste0(dirname,'/','plotssamples_predictAD'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    subsam <- seq(1,dim(distsFA[[acov]])[1], length.out=128)
    tplot(x=agrid, y=t(distsAF[[acov]][,,'AD']),
          col=2, lty=1, lwd=1, alpha=0.875, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          xlab=acov, ylab='probability of AD')
    tplot(x=agrid, y=qdistsAF[[acov]][1,,'AD'],
          col=6, lty=1, lwd=4, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          xlab=acov, ylab='probability of AD', add=T)
    abline(h=0.5, lty=2, lwd=1, col=2)
## legend(x=agrid[1], y=ylim[2]*1, legend=c('87.5% uncertainty on the probability'),
##        col=palette()[c(2)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
##                        )
}
dev.off()


## ## Sample of features of future datapoints
## datasamples <- foreach(asample=1:nrow(parmList$q), .combine=rbind, .packages='nimble', .inorder=F)%dorng%{
##     acluster <- rcat(n=1,prob=parmList$q[asample,])
##     sapply(covNames,function(acov){
##         if(acov %in% realCovs){
##             rnorm(n=1,mean=parmList$meanR[asample,acov,acluster],sd=1/sqrt(parmList$tauR[asample,acov,acluster]))
##         }else if(acov %in% integerCovs){
##             rbinom(n=1,prob=parmList$probI[asample,acov,acluster],size=parmList$sizeI[asample,acov,acluster])
##         }else{
##             nimble::rcat(n=1,prob=c(1-parmList$probB[asample,acov,acluster],parmList$probB[asample,acov,acluster]))-1
##         }
##     })
## }

## pADdatasamples <- samplesF(Y=matrix(1,nrow=1,dimnames=list(NULL,maincov)),
##                            X=datasamples[,otherCovs], parmList=parmList, inorder=F)

## pADdata <- rowMeans(pADdatasamples)
## pADdata <- abs(pADdata-0.5)+0.5
## ## > summary(pADdata)
## ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## ##  0.5004  0.5431  0.5784  0.5863  0.6151  0.7269 

## pADdatasamplestest <- samplesF(Y=matrix(1,nrow=1,dimnames=list(NULL,maincov)),
##                            X=testdata[,..otherCovs], parmList=parmList, inorder=F)

## pADdatatest <- rowMeans(pADdatasamplestest)
## pADdatatest <- abs(pADdatatest-0.5)+0.5
## ## > summary(pADdatatest)
## ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## ##  0.5027  0.5469  0.5836  0.5900  0.6109  0.7221 

## psinglefeatures <- sapply(otherCovs,function(acov){
##     rowMeans(samplesF(Y=matrix(1,nrow=1,dimnames=list(NULL,maincov)),
##              X=datasamples[,acov,drop=F], parmList=parmList, inorder=F))})

## sort(colMeans(abs(psinglefeatures-0.5)+0.5), decreasing=T)
## ## >   RAVLT_immediate  AVDEL30MIN_neuro    AVDELTOT_neuro LRHHC_n_long_log_ 
## ##         0.6580801         0.6536124         0.6254173         0.5931097 
## ##    TRABSCOR_neuro   CATANIMSC_neuro    TRAASCOR_neuro            Apoe4_ 
## ##         0.5889707         0.5838083         0.5706103         0.5447360 
## ##          AGE_log_    ANARTERR_neuro       Gender_num_       GDTOTAL_gds 
## ##         0.5421207         0.5364131         0.5274789         0.5243736 

## entropysort <- sort(colMeans(psinglefeatures*log2(1/psinglefeatures)+(1-psinglefeatures)*log2(1/(1-psinglefeatures))), decreasing=F)
##   ## RAVLT_immediate  AVDEL30MIN_neuro    AVDELTOT_neuro    TRABSCOR_neuro 
##   ##       0.8883023         0.8934159         0.9291837         0.9536632 
##   ##  TRAASCOR_neuro LRHHC_n_long_log_   CATANIMSC_neuro          AGE_log_ 
##   ##       0.9636331         0.9658517         0.9660364         0.9915779 
##   ##          Apoe4_    ANARTERR_neuro       Gender_num_       GDTOTAL_gds 
##   ##       0.9920768         0.9939879         0.9951186         0.9978416 
## message(paste0(names(entropysort),collapse='\n'))




#########################################################
## long-run mutual infos
#########################################################


Xlist <- c(
    as.list(covNames),
    list(otherCovs),
    lapply(otherCovs,function(x){setdiff(otherCovs,x)})
)
names(Xlist) <- c(covNames, 'all', paste0('all_minus_',otherCovs))

## allMI <- samplesMI(Y=maincov, X=Xlist, parmList=parmList, inorder=F, nperf=2^15)
## saveRDS(allMI, paste0(dirname,'/allMI2.rds'))

allMI <- readRDS(paste0(dirname,'/allMI2.rds'))

tquant <- function(xx){yy <- quantile(xx, c(1,4,7)/8, na.rm=T, type=8)
    names(yy) <- c('O1','median','O3')
    yy
}

maincovMI <- allMI[maincov,]

mutualinfo <- tquant(allMI['all',])

condH <- tquant(maincovMI-allMI['all',])

singleMI <- t(sapply(otherCovs, function(acov){
    tquant(allMI[acov,])
}))


singleCH <- t(sapply(otherCovs, function(acov){
    tquant(maincovMI-allMI[acov,])
}))

dropMI <- t(sapply(otherCovs, function(acov){
    tquant(allMI[paste0('all_minus_',acov),])
}))

dropCH <- t(sapply(otherCovs, function(acov){
    tquant(maincovMI-allMI[paste0('all_minus_',acov),])
}))

jointMI <- allMI['all',]
dropMIrel <- t(sapply(otherCovs, function(acov){
    aMI <- allMI[paste0('all_minus_',acov),]
    reldiff <- (1 - aMI/jointMI)*100
    tquant(reldiff)
}))

jointCH <- maincovMI-allMI['all',]
dropCHrel <- t(sapply(otherCovs, function(acov){
    aMI <- maincovMI-allMI[paste0('all_minus_',acov),]
    reldiff <- -(1 - aMI/jointCH)*100
    tquant(reldiff)
}))


######################
#### Save to file ####
outfile2 <- paste0(dirname,'/','results-2_15.txt')
printappr <- function(x,decreasing=T){
    print(signif(x[order(x[,'median'],decreasing=decreasing),], 3))
    ## appr <- cbind(
    ##     ## x,
    ##     signif(x[,1], ceiling(log10(x[,1]/x[,2]))+2),
    ##     signif(x[,2],1))
    ## print(appr[order(x[,1],decreasing=decreasing),])
}
##
sink(outfile2)
cat('Conditional entropy/bit between', maincov, 'and all other features:\n')
printappr(rbind(condH),F)
##
cat('\n\nConditional entropy between', maincov, 'and SINGLE features:\n')
printappr(singleCH,F)
##
cat('\n\nConditional entropy between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
printappr(dropCH,F)
##
cat('\n\nRelative differences between conditional entropy using all features and those using all features minus one, in %:\n')
## relative difference in mutual information without the variate
## printappr(
##     cbind(1-dropMI[,1]/mutualinfo[1],
##           abs(-dropMI[,2]/mutualinfo[1] +
##           mutualinfo[2]*dropMI[,1]/(mutualinfo[1]*mutualinfo[1]))
## #          (dropMI[,2]/dropMI[,1] + mutualinfo[2]/mutualinfo[1]) * (1-dropMI[,1]/mutualinfo[1])
##           )*100
## )
## ##
## cat('\n\n')
print(
    signif(dropCHrel[order(dropCHrel[,2],decreasing=T),],2)
)
##
##
cat('Mutual information/bit between', maincov, 'and all other features:\n')
printappr(rbind(mutualinfo))
##
cat('\n\nMutual information between', maincov, 'and SINGLE features:\n')
printappr(singleMI)
##
cat('\n\nMutual information between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
printappr(dropMI)
##
cat('\n\nRelative differences between mutual information using all features and those using all features minus one, in %:\n')
## relative difference in mutual information without the variate
## printappr(
##     cbind(1-dropMI[,1]/mutualinfo[1],
##           abs(-dropMI[,2]/mutualinfo[1] +
##           mutualinfo[2]*dropMI[,1]/(mutualinfo[1]*mutualinfo[1]))
## #          (dropMI[,2]/dropMI[,1] + mutualinfo[2]/mutualinfo[1]) * (1-dropMI[,1]/mutualinfo[1])
##           )*100
## )
## ##
## cat('\n\n')
print(
    signif(dropMIrel[order(dropMIrel[,2],decreasing=T),],2)
)
##
sink()


########################
#### Plots for talk ####

dropcovs <- paste0('all_minus_',otherCovs)
dropcovsall <- c('all', dropcovs)
shnames <- c(
    sapply(covNames,function(acov){gsub('([^_]+)_.*', '\\1', acov)}),
    all='all',
    sapply(dropcovs,function(acov){gsub('all_(minus_)*([^_]+)_*.*', 'all \\\\ \\2', acov)})
)

thbound <- function(x, N=100){(1/2+(1/2)*sqrt(1-(1-x)^(4/3)))*N}
## thbounderr <- function(x, N=100){num <- 1/2+(1/2)*sqrt(1-(1-x)^(4/3)); c(num*N, 2*sqrt(num*(1-num)*N))}

allothercovs <- setdiff(rownames(allMI),maincov)
maxMI <- max(allMI[allothercovs,])
maxmedianMI <- max(apply(allMI[allothercovs,],1,median))
medianallMI <- apply(allMI,1,median)

ordersingle <- order(medianallMI[otherCovs],decreasing=F)
orderdrop <- order(medianallMI[dropcovs],decreasing=T)
orderdropall <- order(medianallMI[dropcovsall],decreasing=T)

## svnames <- sapply(covNames,function(acov){gsub('([^_]+)_.*', '\\1', acov)})
## ordersingle <- order(singleMI[otherCovs,'median'],decreasing=F)
## ##
## orderdrop <- order(dropMI[,'median'],decreasing=F)
## minusnames <- names(Xlist)[grepl('^all',names(Xlist))]
## minusnames <- minusnames[c(1,orderdrop+1)]
## svminusnames <- sapply(minusnames,function(acov){gsub('all_(minus_)*([^_]+)_.*', 'all \\\\ \\2', acov)})

### Plot of MI of single and discarded features
## pdff(paste0(dirname,'/rankMI_single'))
## tplot(x=medianallMI[otherCovs[ordersingle]], type='p',
##       pch='+',yticks=NA, ylab=NA, xlab='mutual info/Sh', col=1,
##       xlim=c(0,maxMI), ylim=c(0,NA))
## axis(side=3, at=pretty(c(0,maxMI),10), labels=paste0(round(thbound(pretty(c(0,maxMI),10))),'%'), tick=TRUE, lty=1, lwd=0, lwd.ticks=1, col.ticks='#bbbbbb80', cex.axis=1.25, gap.axis=0.25, line=0.5)
## mtext("correct prognoses (TP+TN)", side=3, line=3, cex=1.25)
## for(i in 1:length(otherCovs[ordersingle])){
##     acov <- otherCovs[ordersingle][i]
##     text(x=medianallMI[acov], y=i, labels=shnames[acov], adj=c(0.5,-0.75),xpd=NA, col=1,cex=1)
## }
## tplot(x=medianallMI[dropcovsall[orderdropall]], y=0:length(dropcovs),type='p',
##       pch='+',yticks=NA, ylab=NA, xlab='mutual info/Sh', col=1,
##       xlim=c(0,maxMI),add=T)
## for(i in 1:length(dropcovsall[orderdropall])){
##     acov <- dropcovsall[orderdropall][i]
##     text(x=medianallMI[acov], y=i-1, labels=shnames[acov], adj=c(0.5,-0.75),xpd=NA, col=1,cex=1)
## }
## dev.off()


set.seed(149)
choosesam <- sample(1:ncol(allMI),size=64)
maxMI <- max(allMI[allothercovs,choosesam])
maxy <- 11
plotsingle <- TRUE
plotall <- TRUE
plotdrop <- TRUE
plotunc <- TRUE
### Plot of samples of MI of single and discarded features
pdff(paste0(dirname,'/rankMI_single',plotsingle,'_all',plotall,'_drop',plotdrop,'_unc',plotunc))
if(plotsingle){
    if(plotunc){
tplot(x=allMI[otherCovs[ordersingle],choosesam], type='l',
      pch='+',yticks=NA, ylab=NA, xlab='mutual info/Sh', lty=1,col=5,lwd=1,alpha=0.5,
      xlim=c(0,maxMI), ylim=c(0,maxy))
}
tplot(x=medianallMI[otherCovs[ordersingle]], type='p',
      pch=16, cex=1.25,yticks=NA, ylab=NA, xlab='mutual info/Sh', col=1,
      xlim=c(0,maxMI), ylim=c(0,maxy), add=plotunc)
if(FALSE){
    tplot(x=medianallMI[otherCovs[ordersingle]], type='l',
      pch=16, cex=1,yticks=NA, ylab=NA, xlab='mutual info/Sh', col=1, lty=1, lwd=2,
      xlim=c(0,maxMI), ylim=c(0,maxy), add=T)
}
    for(i in 1:length(otherCovs[ordersingle])){
    acov <- otherCovs[ordersingle][i]
    text(x=medianallMI[acov], y=i, labels=shnames[acov], adj=c(0.5,-0.7),xpd=NA, col='#000000',cex=1)
}
}
##
if(plotdrop){
    if(plotunc){
tplot(x=allMI[dropcovsall[orderdropall],choosesam], y=0:length(dropcovs),type='l',
      pch='+',yticks=NA, ylab=NA, xlab='mutual info/Sh', lty=1,col=2,lwd=1,alpha=0.5,
      xlim=c(0,maxMI), ylim=c(0,maxy), add=plotsingle)
    }
tplot(x=medianallMI[dropcovs[orderdrop]], y=1:length(dropcovs), type='p',
      pch=16, cex=1.25,yticks=NA, ylab=NA, xlab='mutual info/Sh', col=6,
      xlim=c(0,maxMI), ylim=c(0,maxy), add=plotsingle)
if(FALSE){
    tplot(x=medianallMI[dropcovs[orderdrop]], y=1:length(dropcovs), type='l',
      pch=16, cex=1,yticks=NA, ylab=NA, xlab='mutual info/Sh', col=6, lty=1, lwd=2,
      xlim=c(0,maxMI), ylim=c(0,maxy), add=T)
    }
for(i in 1:length(dropcovs[orderdrop])){
    acov <- dropcovs[orderdrop][i]
    text(x=medianallMI[acov], y=i, labels=shnames[acov], adj=c(0.5,-0.7),xpd=NA, col='#000000',cex=1)
}
}
if(plotall | plotdrop){
tplot(x=medianallMI['all'], y=0, type='p',
      pch=10, cex=2,yticks=NA, ylab=NA, xlab='mutual info/Sh', col=3,
      xlim=c(0,maxMI), ylim=c(0,maxy), add=(plotdrop|plotsingle))
text(x=medianallMI['all'], y=0, labels='all features', adj=c(0.5,-0.7),xpd=NA, col=3,cex=1.25)
}
##
axis(side=3, at=pretty(c(0,maxMI),10), labels=paste0(round(thbound(pretty(c(0,maxMI),10))),'%'), tick=TRUE, lty=1, lwd=0, lwd.ticks=1, col.ticks='#bbbbbb80', cex.axis=1.25, gap.axis=0.25, line=0.5)
mtext("max achievable correct prognoses (TP+TN)", side=3, line=3, cex=1.25)
dev.off()





#######################################################
#######################################################
#######################################################
#### Old version
#######################################################
Xlist <- c(
    as.list(covNames),
    list(otherCovs),
    lapply(otherCovs,function(x){setdiff(otherCovs,x)})
)
names(Xlist) <- c(covNames, 'all', paste0('all_minus_',otherCovs))

## allMI <- samplesMI(Y=maincov, X=Xlist, parmList=parmList, inorder=F, nperf=2^13)
## saveRDS(allMI, paste0(dirname,'/allMI.rds'))
allMI <- readRDS(paste0(dirname,'/allMI.rds'))

tquant <- function(xx){yy <- quantile(xx, c(1,4,7)/8, na.rm=T, type=8)
    names(yy) <- c('O1','median','O3')
    yy
}

maincovMI <- allMI[maincov,]

mutualinfo <- tquant(allMI['all',])

condH <- tquant(maincovMI-allMI['all',])

singleMI <- t(sapply(otherCovs, function(acov){
    tquant(allMI[acov,])
}))


singleCH <- t(sapply(otherCovs, function(acov){
    tquant(maincovMI-allMI[acov,])
}))

dropMI <- t(sapply(otherCovs, function(acov){
    tquant(allMI[paste0('all_minus_',acov),])
}))

dropCH <- t(sapply(otherCovs, function(acov){
    tquant(maincovMI-allMI[paste0('all_minus_',acov),])
}))

jointMI <- allMI['all',]
dropMIrel <- t(sapply(otherCovs, function(acov){
    aMI <- allMI[paste0('all_minus_',acov),]
    reldiff <- (1 - aMI/jointMI)*100
    tquant(reldiff)
}))

jointCH <- maincovMI-allMI['all',]
dropCHrel <- t(sapply(otherCovs, function(acov){
    aMI <- maincovMI-allMI[paste0('all_minus_',acov),]
    reldiff <- -(1 - aMI/jointCH)*100
    tquant(reldiff)
}))


######################
#### Save to file ####
outfile2 <- paste0(dirname,'/','results.txt')
printappr <- function(x,decreasing=T){
    print(signif(x[order(x[,'median'],decreasing=decreasing),], 3))
    ## appr <- cbind(
    ##     ## x,
    ##     signif(x[,1], ceiling(log10(x[,1]/x[,2]))+2),
    ##     signif(x[,2],1))
    ## print(appr[order(x[,1],decreasing=decreasing),])
}
##
sink(outfile2)
cat('Conditional entropy/bit between', maincov, 'and all other features:\n')
printappr(rbind(condH),F)
##
cat('\n\nConditional entropy between', maincov, 'and SINGLE features:\n')
printappr(singleCH,F)
##
cat('\n\nConditional entropy between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
printappr(dropCH,F)
##
cat('\n\nRelative differences between conditional entropy using all features and those using all features minus one, in %:\n')
## relative difference in mutual information without the variate
## printappr(
##     cbind(1-dropMI[,1]/mutualinfo[1],
##           abs(-dropMI[,2]/mutualinfo[1] +
##           mutualinfo[2]*dropMI[,1]/(mutualinfo[1]*mutualinfo[1]))
## #          (dropMI[,2]/dropMI[,1] + mutualinfo[2]/mutualinfo[1]) * (1-dropMI[,1]/mutualinfo[1])
##           )*100
## )
## ##
## cat('\n\n')
print(
    signif(dropCHrel[order(dropCHrel[,2],decreasing=T),],2)
)
##
##
cat('Mutual information/bit between', maincov, 'and all other features:\n')
printappr(rbind(mutualinfo))
##
cat('\n\nMutual information between', maincov, 'and SINGLE features:\n')
printappr(singleMI)
##
cat('\n\nMutual information between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
printappr(dropMI)
##
cat('\n\nRelative differences between mutual information using all features and those using all features minus one, in %:\n')
## relative difference in mutual information without the variate
## printappr(
##     cbind(1-dropMI[,1]/mutualinfo[1],
##           abs(-dropMI[,2]/mutualinfo[1] +
##           mutualinfo[2]*dropMI[,1]/(mutualinfo[1]*mutualinfo[1]))
## #          (dropMI[,2]/dropMI[,1] + mutualinfo[2]/mutualinfo[1]) * (1-dropMI[,1]/mutualinfo[1])
##           )*100
## )
## ##
## cat('\n\n')
print(
    signif(dropMIrel[order(dropMIrel[,2],decreasing=T),],2)
)
##
sink()


dropMIq <- t(sapply(otherCovs, function(acov){
    aMI <- allMI[paste0('all_minus_',acov),]
    quantile(aMI, c(1/2, c(1,7)/8))
}))
orderdrop <- order(dropMIq[,1],decreasing=T)

histosmi <- apply(allMI,1,function(aMI){thist(aMI,n=16)})
maxsmi <- sapply(histosmi,function(ahis){max(ahis$density)})
##
## tplot(x=list(histosmi[['FAQ']]$mids, histosmi[['GDTOTAL_gds']]$mids),
##       y=list(histosmi[['FAQ']]$density, histosmi[['GDTOTAL_gds']]$density),
##       border=NA)
##
minusnames <- names(Xlist)[grepl('^all',names(Xlist))]
minusnames <- minusnames[c(1,orderdrop+1)]

set.seed(149)
choosesam <- sample(1:ncol(allMI),size=64)
pdff('_justtestsMI')
tplot(x=allMI[minusnames,choosesam],lty=1,col=3,lwd=1,alpha=0.5,xlim=c(NA,NA),ylabels=minusnames)
tplot(x=allMI[minusnames,choosesam],col=3,alpha=0.5,type='p',pch=20,add=T,xgrid=F,ygrid=F,cex=1,lwd=1)
dev.off()

##
maxmiminus <- round(max(maxsmi[minusnames]))
vspace <- 1
seqbases <- (0:length(minusnames))*(maxmiminus+vspace)
names(seqbases) <- minusnames
##
pdff('_justtestsMI')
tplot(x=lapply(minusnames,function(xx){histosmi[[xx]]$mids}),
      y=lapply(minusnames,function(xx){histosmi[[xx]]$density*8+seqbases[xx]}),
      lty=1,lwd=2,xlim=c(0,0.4))
dev.off()

pdff('_justtestsMI')
tplot(y=lapply(histosmi[!grepl('^all',names(histosmi))],function(xx){xx$mids}),
      x=lapply(histosmi[!grepl('^all',names(histosmi))],function(xx){xx$density*2+seqbases}),
      lty=1,lwd=2)
dev.off()








## YX <- samplesX(nperf=16, parmList=parmList)
## attr(YX, 'rng') <- NULL
## saveRDS(YX, paste0(dirname,'/YX.rds'))
## probYX <- rowMeans(samplesF(Y=YX, parmList=parmList),na.rm=T)
## saveRDS(probYX, paste0(dirname,'/probYX.rds'))
## probY <- rowMeans(samplesF(Y=YX[,maincov,drop=F], parmList=parmList), na.rm=T)
## saveRDS(probY, paste0(dirname,'/probY.rds'))
## probX <- rowMeans(samplesF(Y=YX[,otherCovs,drop=F], parmList=parmList), na.rm=T)
## saveRDS(probX, paste0(dirname,'/probX.rds'))
## ##mutualinfo <- mean(log2(probYX/(probY*probX)), na.rm=T)
## mutualinfo <- c(mean(log2(probYX/(probY*probX)), na.rm=T),
##                 sd(log2(probYX/(probY*probX)), na.rm=T)/errcorr)
## ## > mutualinfo
## ## [1] 0.2073521

## ## probD <- rowMeans(samplesF(Y=YX[,c(integerCovs,binaryCovs),drop=F], parmList=parmList), na.rm=T)

## ## entropy <- mean(log2(1/probD))

## singlemi <- t(sapply(otherCovs, function(acov){
##     probsingle <- rowMeans(samplesF(Y=YX[,acov,drop=F], parmList=parmList), na.rm=T)
##     saveRDS(probsingle,paste0(dirname,'/probsingle_',acov,'.rds'))
##     probjoint <- rowMeans(samplesF(Y=YX[,c(maincov,acov),drop=F], parmList=parmList), na.rm=T)
##     saveRDS(probjoint,paste0(dirname,'/probjoint_',acov,'.rds'))
## ##    mean(log2(probjoint/(probY*probsingle)), na.rm=T)
##     c(mean(log2(probjoint/(probY*probsingle)), na.rm=T),
##       sd(log2(probjoint/(probY*probsingle)), na.rm=T)/errcorr)
## }))
## ## > cbind(sort(singlemi,decreasing=T))
## ##                         [,1]
## ## AVDEL30MIN_neuro 0.102515427
## ## RAVLT_immediate  0.094834522
## ## FAQ              0.084784968
## ## AVDELTOT_neuro   0.058414123
## ## TRABSCOR_neuro   0.035608569
## ## LRHHC_n_long_log 0.031327181
## ## CATANIMSC_neuro  0.028329753
## ## TRAASCOR_neuro   0.027374218
## ## LRLV_n_long_log  0.006005771
## ## AGE_log          0.004534418
## ## GDTOTAL_gds      0.003293993
## ## Gender_num_      0.001216298

## ## singlemi2 <- sapply(otherCovs, function(acov){
## ##     probcondacov <- rowMeans(samplesF(Y=xcond[,2,drop=F], X=YX[,acov,drop=F], parmList=parmList), na.rm=T)
## ## mean(probcondacov*log2(probcondacov/mean(probcondacov,na.rm=T)) + (1-probcondacov)*log2((1-probcondacov)/mean(1-probcondacov,na.rm=T)), na.rm=T)
## ## })

## ## sink(outfile,append=T)
## ## cat('Mutual information between', maincov, 'and SINGLE features, second method:\n')
## ## print(cbind(sort(singlemi2,decreasing=T)))
## ## sink()

## dropMI <- t(sapply(otherCovs, function(acov){
##     probjoint <- rowMeans(samplesF(Y=YX[,setdiff(colnames(YX),acov),drop=F], parmList=parmList),na.rm=T)
##     saveRDS(probjoint,paste0(dirname,'/probjointminus_',acov,'.rds'))
##     probsingle <- rowMeans(samplesF(Y=YX[,setdiff(otherCovs,acov),drop=F], parmList=parmList), na.rm=T)
##     saveRDS(probsingle,paste0(dirname,'/probsingleminus_',acov,'.rds'))
## ##    mean(log2(probjoint/(probY*probsingle)), na.rm=T)
##     c(mean(log2(probjoint/(probY*probsingle)), na.rm=T),
##       sd(log2(probjoint/(probY*probsingle)), na.rm=T)/errcorr)

## }))

## cbind(sort(dropMI,decreasing=F))
## >                       [,1]
## FAQ              0.1827049
## TRABSCOR_neuro   0.1885598
## TRAASCOR_neuro   0.1988626
## RAVLT_immediate  0.2049027
## AVDEL30MIN_neuro 0.2050292
## AVDELTOT_neuro   0.2053379
## LRHHC_n_long_log 0.2057986
## CATANIMSC_neuro  0.2061171
## LRLV_n_long_log  0.2065904
## Gender_num_      0.2068082
## GDTOTAL_gds      0.2073682
## AGE_log          0.2074283


## cbind(sort(1-dropMI/mutualinfo,decreasing=T))*100
## FAQ              11.886637654
## TRABSCOR_neuro    9.062990024
## TRAASCOR_neuro    4.094235113
## RAVLT_immediate   1.181263659
## AVDEL30MIN_neuro  1.120273793
## AVDELTOT_neuro    0.971391568
## LRHHC_n_long_log  0.749197462
## CATANIMSC_neuro   0.595605535
## LRLV_n_long_log   0.367346966
## Gender_num_       0.262317809
## GDTOTAL_gds      -0.007732331
## AGE_log          -0.036716825

## outfile2 <- paste0(dirname,'/','results_err.txt')
## printappr <- function(x){
##     appr <- cbind(
##         x,
##         signif(x[,1], ceiling(log10(x[,1]/x[,2]))+2),
##         signif(x[,2],1))
##     print(appr[order(x[,1],decreasing=T),])
## }
## ##
## sink(outfile2)
## cat('Mutual information/bit between', maincov, 'and all other features:\n')
## printappr(rbind(mutualinfo))
## ##
## cat('\n\nMutual information between', maincov, 'and SINGLE features:\n')
## printappr(singlemi)
## ##
## cat('\n\nMutual information between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
## printappr(dropMI)
## ##
## cat('\n\nRelative differences between mutual information using all features and those using all features minus one, in %:\n')
## ## relative difference in mutual information without the variate
## printappr(
##     cbind(1-dropMI[,1]/mutualinfo[1],
##           abs(-dropMI[,2]/mutualinfo[1] +
##           mutualinfo[2]*dropMI[,1]/(mutualinfo[1]*mutualinfo[1]))
## #          (dropMI[,2]/dropMI[,1] + mutualinfo[2]/mutualinfo[1]) * (1-dropMI[,1]/mutualinfo[1])
##           )*100
## )
## ##
## sink()




## sink(outfile)
## cat('Mutual information between', maincov, 'and all other features:\n',
## mutualinfo, 'bit')
## ## 0.222 bit

## cat('\n\nMutual information between', maincov, 'and SINGLE features:\n')
## print(cbind(sort(singlemi,decreasing=T)))

## cat('\n\nMutual information between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
## print(cbind(sort(dropmis,decreasing=F)))
## ##
## cat('\n\nRelative differences between mutual information using all features and those using all features minus one, in %:\n')
## ## relative difference in mutual information without the variate
## print(cbind(sort(1-dropmis/mutualinfo,decreasing=T))*100)
## ##
## sink()

##print(round(cbind(sort(1-dropmis/mutualinfo,decreasing=T))*100))


## #########################################################
## ## Mutual info for next prediction - combined
## #########################################################
## dirname1 <- 'FAQposteriorAc4_-V12-D708-K75-I1024'
## dirname2 <- 'FAQposteriorAc5_-V12-D708-K75-I1024'
## dirname <- '.'

## errcorr <- sqrt(16*128*2)

## YX <- rbind(readRDS(paste0(dirname1,'/YX.rds')), readRDS(paste0(dirname2,'/YX.rds')))
## probYX <- c(readRDS(paste0(dirname1,'/probYX.rds')), readRDS(paste0(dirname2,'/probYX.rds')))
## probY <- c(readRDS(paste0(dirname1,'/probY.rds')), readRDS(paste0(dirname2,'/probY.rds')))
## probX <- c(readRDS(paste0(dirname1,'/probX.rds')), readRDS(paste0(dirname2,'/probX.rds')))
## ##mutualinfo <- mean(log2(probYX/(probY*probX)), na.rm=T)
## sequ <- log2(probYX/(probY*probX))
## mutualinfo <- c(mean(sequ, na.rm=T),
##                 sd(sequ, na.rm=T)/sqrt(LaplacesDemon::ESS(sequ)))


## singlemi <- t(sapply(otherCovs, function(acov){
##     probsingle <- c(readRDS(paste0(dirname1,'/probsingle_',acov,'.rds')), readRDS(paste0(dirname2,'/probsingle_',acov,'.rds')))
##     probjoint <- c(readRDS(paste0(dirname1,'/probjoint_',acov,'.rds')), readRDS(paste0(dirname2,'/probjoint_',acov,'.rds')))
##     ##
##     sequ <- log2(probjoint/(probY*probsingle))
##     ## print(LaplacesDemon::ESS(sequ))
##     c(mean(sequ, na.rm=T),
##       sd(sequ, na.rm=T)/sqrt(LaplacesDemon::ESS(sequ)))
## }))


## dropmis <- t(sapply(otherCovs, function(acov){
##     probjoint <- c(readRDS(paste0(dirname1,'/probjointminus_',acov,'.rds')), readRDS(paste0(dirname2,'/probjointminus_',acov,'.rds')))
##     probsingle <- c(readRDS(paste0(dirname1,'/probsingleminus_',acov,'.rds')), readRDS(paste0(dirname2,'/probsingleminus_',acov,'.rds')))
##     ##
##     sequ <- log2(probjoint/(probY*probsingle))
##     ##print(LaplacesDemon::ESS(sequ))
##     c(mean(sequ, na.rm=T),
##       sd(sequ, na.rm=T)/sqrt(LaplacesDemon::ESS(sequ)),
##       cov(sequ,
##           log2(probYX/(probY*probX)), use='complete.obs')/(LaplacesDemon::ESS(sequ*log2(probYX/(probY*probX))))
##       )
## }))


## outfile2 <- paste0('results_FAQcombined_err2.txt')
## printappr <- function(x){
##     appr <- cbind(
##         x,
##         signif(x[,1], ceiling(log10(x[,1]/x[,2]))+2),
##         signif(x[,2],1))
##     print(appr[order(x[,1],decreasing=T),])
## }
## ##
## sink(outfile2)
## cat('Mutual information/bit between', maincov, 'and all other features:\n')
## printappr(rbind(mutualinfo))
## ##
## cat('\n\nMutual information between', maincov, 'and SINGLE features:\n')
## printappr(singlemi)
## ##
## cat('\n\nMutual information between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
## printappr(dropmis)
## ##
## cat('\n\nRelative differences between mutual information using all features and those using all features minus one, in %:\n')
## ## relative difference in mutual information without the variate
## printappr(
##     cbind(1-dropmis[,1]/mutualinfo[1],
##           abs(-dropmis[,2]/mutualinfo[1] +
##           mutualinfo[2]*dropmis[,1]/(mutualinfo[1]*mutualinfo[1]))
## #          (dropmis[,2]/dropmis[,1] + mutualinfo[2]/mutualinfo[1]) * (1-dropmis[,1]/mutualinfo[1])
##           )*100
## )
## printappr(
##     cbind(1-dropmis[,1]/mutualinfo[1],
##           abs(dropmis[,1]/mutualinfo[1])*sqrt(
##                                        (dropmis[,2]/dropmis[,1])^2 +
##                                        (mutualinfo[2]/mutualinfo[1])^2 -
##                                        2*dropmis[,3]/(dropmis[,1]*mutualinfo[1])
##                                )
##           )*100
## )
## printappr(
##     cbind(1-dropmis[,1]/mutualinfo[1],
##           abs(dropmis[,1]/mutualinfo[1])*sqrt(
##                                        (dropmis[,2]/dropmis[,1])^2 +
##                                        (mutualinfo[2]/mutualinfo[1])^2 
##                                )
##           )*100
## )
## ##
## sink()




## sink(outfile)
## cat('Mutual information between', maincov, 'and all other features:\n',
## mutualinfo, 'bit')
## ## 0.222 bit

## cat('\n\nMutual information between', maincov, 'and SINGLE features:\n')
## print(cbind(sort(singlemi,decreasing=T)))

## cat('\n\nMutual information between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
## print(cbind(sort(dropmis,decreasing=F)))
## ##
## cat('\n\nRelative differences between mutual information using all features and those using all features minus one, in %:\n')
## ## relative difference in mutual information without the variate
## print(cbind(sort(1-dropmis/mutualinfo,decreasing=T))*100)
## ##
## sink()

## ##print(round(cbind(sort(1-dropmis/mutualinfo,decreasing=T))*100))



###########################
## Exploration on test data
###########################
## Load test file
testdata <- fread(testdatafile, sep=',')
testdata <- testdata[,..covNames]

xcond <- matrix(0:1,ncol=2,dimnames=list(NULL,rep(maincov,2)))
diseasenames <- c('MCI', 'AD')
## subgroup=1 is AD
## subgroup=0 is MCI

## predictive probabilities with uncertainty
ipredictions0 <- samplesF(X=rbind(xcond[,1]), Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)
ipredictions1 <- samplesF(X=rbind(xcond[,2]), Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)

## predictive probabilities via Bayes's theorem
bpredictions <- foreach(adatum=1:nrow(ipredictions0), .combine=rbind)%do%{
    dist <- ipredictions1[adatum,]/(ipredictions0[adatum,]+ipredictions1[adatum,])
}


## direct predictive probabilities 
predictions <- samplesF(Y=rbind(xcond[,2]), X=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)
## predictionsj <- samplesF(Y=cbind(rbind(xcond[,2]),as.matrix(testdata[,..otherCovs])), parmList=parmList, inorder=T)
## predictionsx <- samplesF(Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)
## predictions <- predictionsj/predictionsx
##


##
sink(outfile, append=T)

cat('\n\nAverage uncertainty in test-set predictions, direct:\n',
    mean(rowMeans(predictions,na.rm=T)*log2(1/rowMeans(predictions,na.rm=T))),
    'bit\n')
## > [1] 0.428 bit
cat('SD of uncertainty in test-set predictions, direct:\n',
    sd(rowMeans(predictions,na.rm=T)*log2(1/rowMeans(predictions,na.rm=T)),na.rm=T),
    'bit\n')
## > [1] 0.132 bit
## mean(dbernoulli(x=testdata[[maincov]], prob=rowMeans(predictions,na.rm=T), log=TRUE))
##
cat('\n')
cat('Average uncertainty in test-set predictions, via Bayes:\n',
    mean(rowMeans(bpredictions,na.rm=T)*log2(1/rowMeans(bpredictions,na.rm=T)),na.rm=T),
    'bit\n')
## > [1] 0.420 bit
cat('SD of uncertainty in test-set predictions, via Bayes:\n',
    sd(rowMeans(bpredictions,na.rm=T)*log2(1/rowMeans(bpredictions,na.rm=T)),na.rm=T),
    'bit\n')
## > [1] 0.108 bit
##mean(dbernoulli(x=testdata[[maincov]], prob=rowMeans(bpredictions,na.rm=T), log=TRUE))

sink()



confm <- function(preds, thr=0.5){
c( TP=sum(rowMeans(preds,na.rm=T)>=thr & testdata[[maincov]]==1),
    FP=sum(rowMeans(preds,na.rm=T)>=thr & testdata[[maincov]]==0),
    TN=sum(rowMeans(preds,na.rm=T)<thr & testdata[[maincov]]==0),
    FN=sum(rowMeans(preds,na.rm=T)<thr & testdata[[maincov]]==1)
  )
}

sink(outfile,append=T)

cat('\n\nConfusion matrix test-set, threshold 0.5, direct prediction:\n')
confm(predictions)
## TP FP TN FN 
## 41 27 47 24 
cat('\nConfusion matrix test-set, threshold 0.5, via Bayes:\n')
confm(bpredictions)
## TP FP TN FN 
## 44 32 42 21 

sink()

## tplot(x=agrid <- seq(0,1,length.out=128),
##       y=t(sapply(agrid,confm)))
## legend('top',legend=names(confm()),lty=1:4,col=palette(),bty='n')

## tplot(x=agrid <- seq(0,1,length.out=256),
##       y=colSums(sapply(agrid,confm)[c('FP','FN'),]))


## tplot(x=agrid <- seq(0,1,length.out=256),
##       y=colSums(sapply(agrid,confm)[c('TP','TN'),]))


## tplot(x=agrid <- seq(0,1,length.out=128),
##       y=t(sapply(agrid,bconfm)))
## legend('top',legend=names(bconfm()),lty=1:4,col=palette(),bty='n')

## tplot(x=agrid <- seq(0,1,length.out=256),
##       y=colSums(sapply(agrid,bconfm)[c('FP','FN'),]))

## tplot(x=agrid <- seq(0,1,length.out=256),
##       y=colSums(sapply(agrid,bconfm)[c('TP','TN'),]))

## plot of predictive probabilities for test data
pdff(paste0(dirname,'/','predictions_testset2'))
for(adatum in 1:nrow(testdata)){
    truev <- testdata[[maincov]][adatum]+1
    aprob <- predictions[adatum,]
    meanprob <- mean(aprob, na.rm=T)
    tcol <- 2
    if((meanprob>=0.5 && truev==2) || (meanprob<=0.5 && truev==1)){tcol <- 1}
    uncprob <- quant(aprob, c(1,15)/16)
    histo <- thist(aprob)
    tplot(x=histo$breaks, y=histo$density,
          xlim=c(0,1), ylim=c(0,NA), col=3, 
          xlab='probability of AD', ylab='density')
    abline(v=meanprob, lty=1, lwd=3, col=3)
    legend(x=0.25,y=max(histo$density)*1.15,legend=c(
                         paste0('probability of AD between [',signif(uncprob[1],2),', ',
                                signif(uncprob[2],2),']')),
           cex=1.5, bty='n', xpd=T)
    ## legend('topleft', legend=c(
    ##                       paste0('probability of AD between\n[',signif(uncprob[1],2),', ',signif(uncprob[2],2),']')
    ##                       ),
    ##        cex=1.5, bty='n')
    legend(if(truev==2){pleg <- 'right'}else{pleg <- 'left'},
           legend=c( paste0('true outcome:\n',diseasenames[truev])),
           col=truev, cex=1.5, bty='n')
}
dev.off()

## plot of predictive probabilities for test data using Bayes's theorem
## pdff('predictions_bayes_testset2')
## for(adatum in 1:nrow(testdata)){
##     truev <- testdata[[maincov]][adatum]+1
##     aprob <- bpredictions[adatum,]
##     meanprob <- mean(aprob)
##     tcol <- 2
##     if((meanprob>=0.5 && truev==2) || (meanprob<=0.5 && truev==1)){tcol <- 1}
##     uncprob <- quant(aprob, c(1,15)/16)
##     histo <- thist(aprob)
##     tplot(x=histo$breaks, y=histo$density,
##           xlim=c(0,1), ylim=c(0,NA), col=7, 
##           xlab='uncertainty over the probability of AD', ylab='density')
##     abline(v=mean(aprob), lty=1, lwd=3, col=tcol)
##     abline(v=0.5, lty=3, lwd=2, col=3)
##     legend('topleft', legend=c(
##                           paste0('true outcome: ',diseasenames[truev]),
##                           paste0('probability of AD: ',signif(meanprob*100,2),'%'),
##                           paste0('87.5% uncertainty:\n[',signif(uncprob[1],2),', ',signif(uncprob[2],2),']')
##                           ),
##            cex=1.5, bty='n')
## }
## dev.off()

## Comparison of direct and indirect predictive probabilities
pdff(paste0(dirname,'/','predictions_testset_compbayes2'))
for(adatum in 1:nrow(testdata)){
    truev <- testdata[[maincov]][adatum]+1
    aprob1 <- predictions[adatum,]
    meanprob1 <- mean(aprob, na.rm=T)
    uncprob1 <- quant(aprob, c(1,15)/16)
    histo1 <- thist(aprob)
    aprob <- bpredictions[adatum,]
    meanprob2 <- mean(aprob, na.rm=T)
    uncprob2 <- quant(aprob, c(1,15)/16)
    histo2 <- thist(aprob)
    ymax <- max(histo1$density,histo2$density)
    tplot(x=histo1$breaks, y=histo1$density,
          xlim=c(0,1), ylim=c(0,ymax), col=3, 
          xlab='probability of AD', ylab='density')
    abline(v=meanprob1, lty=1, lwd=3, col=3)
    tplot(x=histo2$breaks, y=histo2$density,
          xlim=c(0,1), ylim=c(0,NA), col=4, 
          xlab='uncertainty over the probability of AD', ylab='density', add=T)
    abline(v=meanprob2, lty=1, lwd=3, col=4)
    ## legend('topleft', legend=c(
    ##                       paste0('probability of AD between\n[',signif(uncprob[1],2),', ',signif(uncprob[2],2),']')
    ##                       ),
    ##        cex=1.5, bty='n')
    legend(x=0.25,y=ymax*1.2,legend=c(
                         paste0('direct prediction [',signif(uncprob1[1],2),', ',
                                signif(uncprob1[2],2),']'),
                         paste0('prediction via bayes [',signif(uncprob2[1],2),', ',
                                signif(uncprob2[2],2),']')),
           col=c(3,4), lty=1, lwd=3, cex=1.5, bty='n', xpd=T)
    legend(if(truev==2){pleg <- 'right'}else{pleg <- 'left'},
           legend=c( paste0('true outcome:\n',diseasenames[truev])),
           col=truev, cex=1.5, bty='n')
}
dev.off()



## #########################################################
## ## 2D plots
## #########################################################

## buildgrid <- function(X, lgrid=128){
##         if(X %in% realCovs){
##             rgx <- range(alldata[[X]])
##             rgx <- rgx + c(-1,1) * diff(rgx)/4
##             lgrid <- seq(rgx[1], rgx[2], length.out=lgrid)
##         }else if(X %in% integerCovs){
##             rgx <- range(alldata[[X]])
##             rgx <- round(rgx + c(-1,1) * diff(rgx)/4)
##             rgx[1] <- max(rgx[1], variateinfo[variate==X,min])
##             rgx[2] <- min(rgx[2], variateinfo[variate==X,max])
##             if(diff(rgx)<lgrid){lgrid <- rgx[1]:rgx[2]}
##             else{lgrid <- round(seq(rgx[1], rgx[2], length.out=lgrid))}
##         }else{
##             lgrid <- 0:1
##         }
##         lgrid
## }

## acov2 <- 'AVDEL30MIN_neuro'#integerCovs[1]
## acov1 <- maincov
## ##
## xgrid <- buildgrid(acov1, 128)
## ygrid <- buildgrid(acov2, 128)
## ##
## grid2d <- cbind(rep(xgrid,length(ygrid)), rep(ygrid, each=length(xgrid)))
## colnames(grid2d) <- c(acov1, acov2)
## ##
## nfsamples <- 32
## fsamples2d <- samplesF(Y=grid2d, parmList=parmList, nfsamples=nfsamples, inorder=F)
## ##
## ##dim(fsamples2d) <- c(length(xgrid), length(ygrid), nfsamples)
## ##
## asample <- 1
## ## ax <- min(diff(xgrid)[1], diff(ygrid)[1])/2
## ## ay <- min(diff(xgrid)[1], diff(ygrid)[1])/2
## ax <- diff(xgrid)[1]/2
## if(acov1 %in% realCovs){ xticks <- NULL }else{ xticks <- xgrid }
## ay <- diff(ygrid)[1]/2
## if(acov2 %in% realCovs){ yticks <- NULL }else{ yticks <- ygrid }
## pmax <- max(fsamples2d[,asample])
## ##
## pdff('_test2dplot')
## plot2dF(xygrid=grid2d, fsamples=rowMeans(fsamples2d,na.rm=T))
## #plot2dF(xygrid=grid2d, fsamples=apply(fsamples2d,1,function(x)diff(quant(x,c(1,15)/16,na.rm=T))))
## plot2dF(xygrid=grid2d, fsamples=apply(fsamples2d,1,function(x)IQR(x,type=8,na.rm=T)))
## dev.off()




## pdff('_test2dplot')
## par(mfrow=c(4,8))
## for(asample in 1:nfsamples){
## plot2dF(xygrid=grid2d, fsamples=fsamples2d[,asample], ticks=F, labs=F, mar=c(0,0,0,0))
##     }
## dev.off()



## pdff('_test2dplot')
## plot2dF(xygrid=grid2d, fsamples=fsamples2d[,1])
## dev.off()

## tplot(x=NA, y=NA, xlim=extendrange(xgrid), ylim=extendrange(ygrid), xlab=acov1, ylab=acov2, xticks=xticks, yticks=yticks)
## for(i in 1:nrow(grid2d)){
##     rat <- fsamples2d[i,asample]/pmax
##     polygon(x=grid2d[i,1]+c(-1,1,1,-1)*ax,
##             y=grid2d[i,2]+c(-1,-1,1,1)*ay,
##             border=gray(1-rat), col=gray(1-rat))
## }
## #tplot(x=NA, y=NA, xlim=extendrange(xgrid), ylim=extendrange(ygrid), xlab=acov1, ylab=acov2,add=T)
## dev.off()







## ##
## pdff('_test2dplot')
## tplot(x=NA, y=NA, xlim=extendrange(xgrid), ylim=extendrange(ygrid), xlab=acov1, ylab=acov2)
## for(i in 1:nrow(grid2d)){
##     rat <- sqrt(fsamples2d[i,asample]/pmax)
##     polygon(x=grid2d[i,1]+c(-1,1,1,-1)*ax*rat,
##             y=grid2d[i,2]+c(-1,-1,1,1)*ay*rat,
##             border='white', col='black')
## }
## dev.off()


## plot2DsamplesF <- function(X, Y, parmList, xgrid=128, ygrid=128){
##     if(length(xgrid)==1){
##         if(X %in% realCovs){
##             rgx <- range(alldata[[X]])
##             rgx <- rgx + c(-1,1) * diff(rgx)/4
##             xgrid <- seq(rgx[1], rgx[2], length.out=xgrid)
##         }else if(X %in% integerCovs){
##             rgx <- range(alldata[[X]])
##             rgx <- round(rgx + c(-1,1) * diff(rgx)/4)
##             rgx[1] <- max(rgx[1], thminicovs[X])
##             rgx[2] <- min(rgx[2], thmaxicovs[X])
##             if(diff(rgx)<xgrid){xgrid <- rgx[1]:rgx[2]}
##             else{xgrid <- round(seq(rgx[1], rgx[2], length.out=xgrid))}
##         }else{
##             xgrid <- 0:1
##         }
##         ##
##         if(Y %in% realCovs){
##             rgy <- range(alldata[[Y]])
##             rgy <- rgy + c(-1,1) * diff(rgy)/4
##             ygrid <- seq(rgy[1], rgy[2], length.out=ygrid)
##         }else if(Y %in% integerCovs){
##             rgy <- range(alldata[[Y]])
##             rgy <- round(rgy + c(-1,1) * diff(rgy)/4)
##             rgy[1] <- max(rgy[1], thminicovs[Y])
##             rgy[2] <- min(rgy[2], thmaxicovs[Y])
##             if(diff(rgy)<ygrid){ygrid <- rgy[1]:rgy[2]}
##             else{ygrid <- round(seq(rgy[1], rgy[2], length.out=ygrid))}
##         }else{
##             ygrid <- 0:1
##         }
##         ##
        


        
##         if(length(ygrid)==1){
##             rgx <- range(alldata[[X]])
##             rgx <- rgx + c(-1,1) * diff(rgx)/4
##             xgrid <- seq(rgx[1], rgx[2], length.out=xgrid)
##         }
## }
