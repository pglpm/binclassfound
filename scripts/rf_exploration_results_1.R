## Author: PGL  Porta Mana
## Created: 2022-03-17T14:21:57+0100
## Last-Updated: 2022-04-28T18:29:54+0200
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
library('LaplacesDemon')
## library('extraDistr')
## library('mvtnorm')
## options(bitmapType='cairo')
## pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
## pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
## library('nimble')
#### End custom setup ####
#### FUNCTION TO CALCULATE MUTUAL INFO FROM JOINT DISTRIBUTION
## freqs[S,B] = freq spike count B and stimulus S (one ROW per stimulus)
## The function calculates the conditional frequencies of B|S
## and if requested constructs
## a new joint distribution with equal marginals for S
## Note: don't need to normalize input to mutualinfo
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
accmi <- function(mi){0.5 + 0.5 * sqrt(1 - (1 - mi)^(4/3))}
entropy <- function(freqs, base=2L){
    freqs <- cbind(freqs)
    freqs <- t(t(freqs)/colSums(freqs, na.rm=T))
    c(entropy=colSums(freqs*log2(1/freqs), na.rm=T)/log2(base))
}
condpRows <- function(x){
    t(t(x)/colSums(x))
}
condpCols <- function(x){
    x/rowSums(x)
}

dt  <- as.data.table(read.csv('modCHEMBL205_predictions_CNN.csv',header=TRUE,sep=','))





nsamples <- as.integer(40000)
tausamples <- -log(rgamma(n=nsamples, shape=1/4, rate=(1/0.3)^2))/2
qt <- quantile(tausamples, (1:3)/4)
##
his <- thist(tausamples)
tplot(x=his$mids, y=his$density, xlab='ld-SD', ylab='density',ylim=c(0,NA))
abline(v=0,col=2)
for(i in qt){abline(v=i,col=4,lty=2)}
exp(min(tausamples))























##
class0 <- dt[class==0]
class1 <- dt[class==1]

h0 <- thist(class0$softmax_output0)
h1 <- thist(class1$softmax_output0)
##
hb0 <- thist(class0$sigmoid_output0)
hb1 <- thist(class1$sigmoid_output0)

nbin <- 16+1
##
h0 <- thist(class0$softmax_output0, n=seq(0,1,length.out=nbin))
h1 <- thist(class1$softmax_output0, n=seq(0,1,length.out=nbin))
##
hb0 <- thist(class0$sigmoid_output0, n=seq(0,1,length.out=nbin))
hb1 <- thist(class1$sigmoid_output0, n=seq(0,1,length.out=nbin))
##
## tplot(h0$mids, cbind(h1$density-h0$density, hb1$density-hb0$density))
##
pdff('softmax_vs_prob')
tplot(h0$mids, h0$counts/(h0$counts+h1$counts), xlab='softmax output 0', ylab='probability of class 0', xlim=c(0,1), ylim=c(0,1))
dev.off()

##
jp <- rbind(h0$counts,h1$counts)
jpb <- rbind(hb0$counts,hb1$counts)
##
mutualinfo(jp)
mutualinfo(jpb)

accmi(mutualinfo(jp))
accmi(mutualinfo(jpb))

entropy(rowSums(jp)) - mutualinfo(jp)
entropy(rowSums(jpb)) - mutualinfo(jpb)
##entropy(colSums(jp)) - mutualinfo(jpb)


############################################################
## Write data to file for nonparametric analysis

output0 <- qlogis(dt$sigmoid_output0)
output1 <- qlogis(dt$sigmoid_output1)

outdt <- data.table(item=as.integer(1:nrow(dt)), class=as.integer(dt$class), output0=output0, output1=output1)

fwrite(outdt, 'softmaxdata_test2.csv', sep=',')
##
outdtshuffled <- outdt[sample(1:nrow(outdt)),]
fwrite(outdtshuffled, 'softmaxdata_test2_shuffled.csv', sep=',')



softmax2 <- ((dt$softmax_output0-0.5)*(1-min(dt$softmax_output0)))+0.5
outdt <- data.table(item=as.integer(1:nrow(dt)), class=as.integer(dt$class), logitsoftmax=qlogis(softmax2))
outdt <- outdt[sample(1:nrow(outdt)),]
fwrite(outdt, 'softmaxdata_test.csv', sep=',')

outdtprobs <- samplesF(Y=cbind(class=0), X=cbind(logitsoftmax=qlogis(softmax2)), parmList=parmList, inorder=T)

newdata <- data.table(class=dt$class,
                      softmax_output0=softmax2,
                      probability0=rowMeans(outdtprobs))

fwrite(newdata, 'CHEMBL205_cl_softmax0_prob0.csv', sep=',')






## We consider two binary {0,1} variables X and Y and a population with a given conditional frequency of Y given X

pXgY0 <- rbind(1-0.4375, 0.4375)
pXgY1 <- rbind(0.25, 1-0.25)
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
##
##
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
##
print('distr:')
pYX
print('true:')
pYgX
print('mod:')
pYgX2
##
print('acc true:')
sum(sapply(1:ncol(pYX),function(x)pYX[which.max(pYgX[,x]),x]))
print('acc mod:')
sum(sapply(1:ncol(pYX),function(x)pYX[which.max(pYgX2[,x]),x]))


#########################################################
#########################################################
#########################################################



set.seed(707)
baseversion <- '_presults1'
nclusters <- 64L
niter <- 1024L # iterations AFTER thinning
niter0 <- 1024L*2L
thin <- 1L
nstages <- 0L
ncheckprobs1 <- 16L
ncheckprobs2 <- 8L
maincov <- 'class'
family <- 'Palatino'
ndata <- 8192L #4096L
posterior <- TRUE
##
## stagestart <- 0L # last saved + 1
##
saveinfofile <- 'variate_info2.csv'
datafile <- 'softmaxdata_test2_shuffled.csv'
#64K, 8192D, 1024I: 0.5 h
odata <- fread(datafile, sep=',')
alldata <- odata[1:ndata, ..covNames]
source('functions_mcmc.R')


## baseversion <- paste0(baseversion,'_',mcmcseed,'_')
variateinfo <- fread(saveinfofile, sep=',')
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
restcovs <- setdiff(covNames, maincov)

dirname <- '_presults1-V3-D8192-K64-I1024'
nfromeach <- 2048/128
mcsamples <- foreach(i=1:128, .combine=rbind)%dopar%{
    temp <- readRDS(paste0(dirname, '/_mcsamples-R_presults1_',i,'_0-V3-D8192-K64-I2048.rds'))
    temp[nrow(temp)+1-(nfromeach:1),]
}
parmlist <- mcsamples2parmlist(mcsamples)

parmlistb <- mcsamples2parmlist(rbind(
    readRDS(paste0('_results5_111_-V3-D8192-K64-I1024', '/_mcsamples-R_results5_111_2-V3-D8192-K64-I1024.rds')),
    readRDS(paste0('_results5_222_-V3-D8192-K64-I1024', '/_mcsamples-R_results5_222_3-V3-D8192-K64-I1024.rds'))
    ))

condtrace1 <- logsumsamplesF(Y=as.matrix(alldata[,..maincov]), X=as.matrix(alldata[,..restcovs]), parmList=parmlist, inorder=F)

condtrace2 <- logsumsamplesF(Y=as.matrix(alldata[,..restcovs]), X=as.matrix(alldata[,..maincov]), parmList=parmlist, inorder=F)

condtrace1b <- logsumsamplesF(Y=as.matrix(alldata[,..maincov]), X=as.matrix(alldata[,..restcovs]), parmList=parmlistb, inorder=F)

condtrace2b <- logsumsamplesF(Y=as.matrix(alldata[,..restcovs]), X=as.matrix(alldata[,..maincov]), parmList=parmlistb, inorder=F)


tplot(y=list(condtrace2,condtrace2b))



xgrid <- seq(0, 1, length.out=256)
ygrid <- X2Y[['prediction_lnodds']](xgrid)
##
vpoints <- cbind(prediction_lnodds=ygrid)

pgrid <- samplesF(Y=cbind(class=0), X=vpoints, parmList=parmList, inorder=F)

qgrid <- apply(pgrid,1,function(x){quantile(x, c(1,7)/8)})
##
tplot(x=xgrid, y=rowMeans(pgrid), xlab='RF % output', ylab='probability of class 0', ylim=c(0,1))
polygon(x=c(xgrid,rev(xgrid)), y=c(qgrid[1,],rev(qgrid[2,])), col=paste0(palette()[1],'40'), border=NA)
legend('topleft', legend=c(
                       paste0(paste0(rownames(qgrid),collapse='\u2013'), ' uncertainty')
                   ),
       lty=c('solid'), lwd=c(10),
       col=paste0(palette()[1],c('40')),
       bty='n', cex=1.25)



pgrid0 <- samplesF(X=cbind(class=0), Y=vpoints, parmList=parmList, inorder=F)*Xjacobian[['prediction_lnodds']](xgrid)
##
qgrid0 <- apply(pgrid0,1,function(x){quantile(x, c(1,7)/8)})
##
pgrid1 <- samplesF(X=cbind(class=1), Y=vpoints, parmList=parmList, inorder=F)*Xjacobian[['prediction_lnodds']](xgrid)
##
qgrid1 <- apply(pgrid1,1,function(x){quantile(x, c(1,7)/8)})
##

tplot(x=xgrid, y=cbind(rowMeans(pgrid0),rowMeans(pgrid1)), xlab='RF % output', ylab='probability of output', ylim=c(0,25))
polygon(x=c(xgrid,rev(xgrid)), y=c(qgrid0[1,],rev(qgrid0[2,])), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(qgrid1[1,],rev(qgrid1[2,])), col=paste0(palette()[2],'40'), border=NA)
legend('topleft', legend=c(
                      'Conditional on class 0',
                      'Conditional on class 1',
                       paste0(paste0(rownames(qgrid),collapse='\u2013'), ' uncertainty')
                   ),
       lty=c(1,2,1), lwd=c(3,3,10),
       col=paste0(palette()[c(1,2,7)],c('','','C0')),
       bty='n', cex=1.25)












orange <- c(-1,1)*ceiling(max(abs(as.matrix(odata[,..realCovs]))))

cseq <- seq(orange[1], orange[2], length.out=128)


vpoints <- cbind(output0=rep(cseq, length(cseq)), output1=rep(cseq, each=length(cseq)))


pgrid <- samplesF(Y=cbind(class=1), X=vpoints, parmList=parmlist, inorder=F)

mpgrid <- rowMeans(pgrid)
dim(mpgrid) <- rep(length(cseq), 2)
##
image(z=mpgrid, x=cseq, y=cseq, zlim=c(0,1), col=gray.colors(128,start=1,end=0), xlab=realCovs[1], ylab=realCovs[2])
grid(lty=1,nx=8,ny=8)

spgrid <- apply(pgrid,1,IQR)
dim(spgrid) <- rep(length(cseq), 2)
##
image(z=spgrid, x=cseq, y=cseq, zlim=c(0,1), col=gray.colors(128,start=1,end=0), xlab=realCovs[1], ylab=realCovs[2])
grid(lty=1,nx=8,ny=8)




cnndata <- fread('modCHEMBL205_predictions_CNN.csv', sep=',')
rfdata <- fread('modCHEMBL205_predictions_RF.csv', sep=',')
svmdata <- fread('modCHEMBL205_predictions_SVM.csv', sep=',')

rfdata <- fread('modCHEMBL205_predictions_RF.csv', sep=',')

rfdata$prediction_int <- as.integer(round(rfdata$prediction * 200))
max(abs(rfdata$prediction_int-rfdata$prediction*200))

fwrite(rfdata, 'modCHEMBL205_predictions_RF.csv', sep=',')


rfdata <- fread('modCHEMBL205_predictions_RF.csv', sep=',')

shrink <- function(x){(x-0.5)*(1-2^-10)+0.5}

rfdata$prediction_lnodds <- log(shrink(rfdata$prediction)/(1-shrink(rfdata$prediction)))

fwrite(rfdata, 'modCHEMBL205_predictions_RF.csv', sep=',')

