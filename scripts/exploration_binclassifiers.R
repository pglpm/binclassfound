## Author: PGL  Porta Mana
## Created: 2022-03-17T14:21:57+0100
## Last-Updated: 2022-04-09T13:56:46+0200
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


dt  <- as.data.table(read.csv('CHEMBL205_cl_sigmoid_and_softmax.csv',header=TRUE,sep=','))
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
softmax2 <- ((dt$softmax_output0-0.5)*(1-min(dt$softmax_output0)))+0.5
outdt <- data.table(item=as.integer(1:nrow(dt)), class=as.integer(dt$class), logitsoftmax=qlogis(softmax2))
outdt <- outdt[sample(1:nrow(outdt)),]
fwrite(outdt, 'softmaxdata_test.csv', sep=',')









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
