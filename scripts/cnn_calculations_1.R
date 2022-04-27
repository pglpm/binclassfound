## Author: PGL  Porta Mana
## Created: 2022-03-17T14:21:57+0100
## Last-Updated: 2022-04-27T16:44:05+0200
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
##library('LaplacesDemon')
## library('extraDistr')
## library('mvtnorm')
## options(bitmapType='cairo')
## pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
## pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
## library('nimble')
#### End custom setup ####

source('functions_mcmc.R')
set.seed(707)
baseversion <- '_presults2'
maincov <- 'class'
family <- 'Palatino'
saveinfofile <- 'cnn_variateinfo.csv'
datafile <- 'modCHEMBL205_predictions_CNN.csv'
##
variateinfo <- fread(saveinfofile, sep=',')
covNames <- variateinfo$variate
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- covNames
odata <- fread(datafile, sep=',')
##
realCovs <- covNames[covTypes=='double']
integerCovs <- covNames[covTypes=='integer']
binaryCovs <- covNames[covTypes=='binary']
covNames <- c(realCovs, integerCovs, binaryCovs)
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(odata)}
alldata <- odata[1:ndata, ..covNames]

dirname <- '_presults2-V3-D3588-K64-I2048'

#### Compare diagnostic traces from parallel runs
parrange <- 3:16
snip <- 256:1
traces <- foreach(i=parrange, .combine=rbind)%dopar%{
    temp <- readRDS(paste0(dirname,'/_probtraces-R_presults2_',i,'_1-V3-D3588-K64-I2048.rds'))
    temp[nrow(temp)+1-snip,]*10/log(10)
}
## tplot(y=traces[,3],lty=1)

temp <- traces
dim(temp) <- c(length(snip), length(parrange), 3)
dimnames(temp) <- list(NULL, NULL, colnames(traces))
## tplot(y=temp[,,1],lty=1,lwd=0.5,col=7)
rm(temp)


#### Collect MCMC data
parrange <- 1:16
snip <- 256:1
mcsamples <- foreach(i=parrange, .combine=rbind)%dopar%{
    temp <- readRDS(paste0(dirname,'/_mcsamples-R_presults2_',i,'_1-V3-D3588-K64-I2048.rds'))
    temp[nrow(temp)+1-snip,]
}
parmlist <- mcsamples2parmlist(mcsamples)

oneparmlist <- list(
    q=rbind(c(parmlist$q))/nrow(parmlist$q),
    meanR=array(aperm(parmlist$meanR, c(2,1,3)),
                dim=c(1, dim(parmlist$meanR)[2], length(parmlist$q)),
                dimnames=dimnames(parmlist$meanR)),
    tauR=array(aperm(parmlist$tauR, c(2,1,3)),
                dim=c(1, dim(parmlist$tauR)[2], length(parmlist$q)),
                dimnames=dimnames(parmlist$tauR)),
    probI=array(aperm(parmlist$probI, c(2,1,3)),
                dim=c(1, dim(parmlist$probI)[2], length(parmlist$q)),
                dimnames=dimnames(parmlist$probI)),
    sizeI=array(aperm(parmlist$sizeI, c(2,1,3)),
                dim=c(1, dim(parmlist$sizeI)[2], length(parmlist$q)),
                dimnames=dimnames(parmlist$sizeI)),
    probB=array(aperm(parmlist$probB, c(2,1,3)),
                dim=c(1, dim(parmlist$probB)[2], length(parmlist$q)),
                dimnames=dimnames(parmlist$probB))
)

qlist <- oneparmlist$q
qorder <- order(c(oneparmlist$q), decreasing=TRUE)
qcumsum <- cumsum(qlist[qorder])

qorder <- order(c(oneparmlist$q), decreasing=TRUE)
##
maxclusters <- 2^10
##
fq <- oneparmlist$q[qorder[1:maxclusters]]
fq <- fq/sum(fq)
fmeanR <- oneparmlist$meanR[1,,qorder[1:maxclusters]]
ftauR <- oneparmlist$tauR[1,,qorder[1:maxclusters]]
fprobB <- oneparmlist$probB[1,,qorder[1:maxclusters]]



##########################################################
##########################################################
##########################################################





dt  <- as.data.table(read.csv('modCHEMBL205_predictions_CNN.csv',header=TRUE,sep=','))





nsamples <- as.integer(2^16)
tausamples <- -log10(rgamma(n=nsamples, shape=1/4, rate=1/10^2))/2
qt <- quantile(tausamples, (1:3)/4)
##
his <- thist(tausamples)
tplot(x=his$mids, y=his$density, xlab='ld-SD', ylab='density',ylim=c(0,NA))
abline(v=0,col=2)
for(i in qt){abline(v=i,col=4,lty=2)}
























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

