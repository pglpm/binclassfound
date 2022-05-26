## Author: PGL  Porta Mana
## Created: 2022-03-17T14:21:57+0100
## Last-Updated: 2022-05-26T15:46:24+0200
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
library('future.apply')
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
library('plotly')
#### End custom setup ####

set.seed(707)
baseversion <- '_cnn_new'
maincov <- 'class'
outputcov <- c('output0', 'output1')
family <- 'Palatino'
saveinfofile <- 'cnn_variateinfo.csv'
calibfile <- 'modCHEMBL205_predictions_CNN_test1_calibration.csv'
demofile <- 'modCHEMBL205_predictions_CNN_test2_demonstration.csv'

##
variateinfo <- fread(saveinfofile, sep=',')
covNames <- variateinfo$variate
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- covNames
cdata <- fread(calibfile, sep=',')
##
realCovs <- covNames[covTypes=='double']
integerCovs <- covNames[covTypes=='integer']
binaryCovs <- covNames[covTypes=='binary']
covNames <- c(realCovs, integerCovs, binaryCovs)
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(cdata)}
alldata <- cdata[1:ndata, ..covNames]
##
X2Y <- list()
Xjacobian <- list()
Xrange <- list()
for(avar in realCovs){
    if(is.null(X2Y[[avar]])){X2Y[[avar]] <- identity}
    if(is.null(Xjacobian[[avar]])){Xjacobian[[avar]] <- function(x){1}}
    if(is.null(Xrange[[avar]])){
        datum <- alldata[1:ndata][[avar]]
        Xrange[[avar]] <- range(datum, na.rm=T)+c(-1,1)*IQR(datum, type=8, na.rm=T)
    }
}
##
source('functions_mcmc.R')
dirname <- '_cnn_new-V3-D3589-K64-I2048'

npar <- 16
ntotal <- 1024*4
nskip <- 8
parmlist <- mcsamples2parmlist(
    foreach(i=1:npar, .combine=rbind)%dopar%{
        temp <- readRDS(paste0(dirname, '/_mcsamples-R_cnn_new_',i,'_1-V3-D3589-K64-I2048.rds'))
        if(any(is.na(nrow(temp)+1-rev(seq(1,nrow(temp),by=nskip)[1:(ntotal/npar)])))){print('WARNING! not enough points')}
        temp[nrow(temp)+1-rev(seq(1,nrow(temp),by=nskip)[1:(ntotal/npar)]),]
}
)

#########################################################
## transducer curve p(c | y)
#########################################################

xr <- round(max(abs(range(unlist(Xrange))))) * c(-1,1)
cseq <- seq(xr[1], xr[2], length.out=129)
##
vpoints <- cbind(rep(cseq, length(cseq)), rep(cseq, each=length(cseq)))
colnames(vpoints) <- realCovs

system.time(opgrid <- samplesF(Y=cbind(class=1), X=vpoints, parmList=parmlist, inorder=F))
##

## system.time(opgrid2 <- samplesFp(Y=cbind(class=1), X=vpoints, batchsize=1024, parmList=parmlist, inorder=F))


mpgrid <- rowMeans(opgrid)
dim(mpgrid) <- rep(length(cseq), 2)

fig <- plot_ly(z=t(mpgrid), x=cseq, y=cseq, cmin=0, cmax=1)
fig <- fig %>% add_surface(colors='Blues')
fig <- fig %>% layout(scene = list(
                          xaxis = list(title = "output 0"),
                          yaxis = list(title = "output 1"),
                          zaxis = list(title = 'P(class 1 | outputs)', range=c(0,1)),
                          camera = list(projection = list(type = 'orthographic'))
                      ), title=NA)
fig
#orca(fig, '../transducer_surface_CNN.pdf', scale=1, width=16.5, height=16.5)

softm <- apply(vpoints, 1, function(x){exp(x[2])/sum(exp(x))})
dim(softm) <- rep(length(cseq), 2)
fig <- plot_ly(z=t(softm), x=cseq, y=cseq, cmin=0, cmax=1)
fig <- fig %>% add_surface(colors='Blues')
fig <- fig %>% layout(scene = list(
                          xaxis = list(title = "output 0"),
                          yaxis = list(title = "output 1"),
                          zaxis = list(title = 'softmax', range=c(0,1)),
                          camera = list(projection = list(type = 'orthographic'))
                      ), title=NA)
fig
#orca(fig, '../transducer_surface_CNN.pdf', scale=1, width=16.5, height=16.5)



diagon <- which(vpoints[,1]==-vpoints[,2])
## softm <- apply(vpoints[diagon,], 1, function(x){exp(x[2])/sum(exp(x))})
xgrid <- vpoints[diagon,2]
pdff('../transducer_curve_diagonal_CNN')
tplot(x=xgrid, y=cbind(rowMeans(opgrid[diagon,]), 1- rowMeans(opgrid[diagon,])), xlab=bquote('output 1' == -'output 0'),
##      ylab=expression(p~group('(',class~output,')')),
      ylab=bquote('P'~group('(','class', '.')~group('|', ' output 1',')')),
      ylim=c(0,1), mar=c(4.5,5.75,1,1), cex.axis=2, cex.lab=2,
       lwd=3, family='Palatino')
legend('topleft', c('class 1', 'class 0'), lty=c(1,2), col=c(1,2), lwd=3, bty='n', cex=2)
##
## polygon(x=c(xgrid,rev(xgrid)), y=c(q1grid[diagon],rev(q2grid[diagon])), col=paste0(palette()[1],'40'), border=NA)
## polygon(x=c(xgrid,rev(xgrid)), y=c(1-q2grid[diagon],rev(1-q1grid[diagon])), col=paste0(palette()[2],'40'), border=NA)
dev.off()

diagon <- which(vpoints[,1]==-vpoints[,2])
## softm <- apply(vpoints[diagon,], 1, function(x){exp(x[2])/sum(exp(x))})
xgrid <- vpoints[diagon,2]
pdff('../transducer_curve_diagonal_CNN_prob')
tplot(x=xgrid, y=rowMeans(opgrid[diagon,]), xlab=bquote('output 1' == -'output 0'),
##      ylab=expression(p~group('(',class~output,')')),
      ylab=bquote('P'~group('(','class 1', '.')~group('|', ' output 1',')')),
      ylim=c(0,1), mar=c(4.5,5.75,1,1), cex.axis=2, cex.lab=2,
       lwd=3, family='Palatino')
#legend('topleft', c('class 1', 'class 0'), lty=c(1,2), col=c(1,2), lwd=3, bty='n', cex=2)
##
polygon(x=c(xgrid,rev(xgrid)), y=c(q1grid[diagon],rev(q2grid[diagon])), col=paste0(palette()[1],'40'), border=NA)
## polygon(x=c(xgrid,rev(xgrid)), y=c(1-q2grid[diagon],rev(1-q1grid[diagon])), col=paste0(palette()[2],'40'), border=NA)
dev.off()

tplot(x=xgrid, y=q1grid[diagon])
tplot(x=xgrid, y=q2grid[diagon],add=T)


q1grid <- apply(opgrid,1,function(x){quantile(x, 1/8)})
##dim(q1grid) <- rep(length(cseq), 2)
q2grid <- apply(opgrid,1,function(x){quantile(x, 7/8)})
##dim(q2grid) <- rep(length(cseq), 2)




####
iqrgrid <- apply(opgrid,1,IQR)
dim(iqrgrid) <- rep(length(cseq), 2)

fig <- plot_ly(z=t(iqrgrid), x=cseq, y=cseq, cmin=0, cmax=1)
fig <- fig %>% add_surface(colors='Reds')
fig <- fig %>% layout(scene = list(xaxis = list(title = "output 0"), yaxis = list(title = "output 1"), zaxis = list(title = "interquartile range around probability of class 1", range=c(0,1))))
fig

####
q1grid <- apply(opgrid,1,function(x){quantile(x, 1/8)})
dim(q1grid) <- rep(length(cseq), 2)
q2grid <- apply(opgrid,1,function(x){quantile(x, 7/8)})
dim(q2grid) <- rep(length(cseq), 2)

fig <- plot_ly(x=cseq, y=cseq, cmin=0, cmax=1)
fig <- fig %>% add_surface(z = t(mpgrid), opacity = 0.5, colorscale='Blues')
fig <- fig %>% add_surface(z = t(q1grid), opacity = 0.5, colorscale='Reds')
fig <- fig %>% add_surface(z = t(q2grid), opacity = 0.5, colorscale='Reds')
fig <- fig %>% layout(scene = list(xaxis = list(title = "output 0"), yaxis = list(title = "output 1"), zaxis = list(title = "probability of class 1", range=c(0,1))))
fig


#########################################################
## Put all parameters into one list
#########################################################
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
##
qorder <- order(c(oneparmlist$q), decreasing=FALSE)
shortparmlist <- list(
    q=oneparmlist$q[,qorder, drop=F]/sum(oneparmlist$q[,qorder]),
    meanR=oneparmlist$meanR[,,qorder, drop=F],
    tauR=oneparmlist$tauR[,,qorder, drop=F],
    probI=oneparmlist$probI[,,qorder, drop=F],
    sizeI=oneparmlist$sizeI[,,qorder, drop=F],
    probB=oneparmlist$probB[,,qorder, drop=F]
)
fwrite(data.table(w=c(shortparmlist$q),
               p=c(shortparmlist$probB),
               mu0=c(shortparmlist$meanR[,'output0',]),
               sigma0=1/sqrt(c(shortparmlist$tauR[,'output0',])),
               mu1=c(shortparmlist$meanR[,'output1',]),
               sigma1=1/sqrt(c(shortparmlist$tauR[,'output1',]))
               ), '_CNN_transducer_parameters.csv', sep=',')



#########################################################
## Calculation of utility yields on demonstration set
#########################################################
ddata <- fread(demofile, sep=',')
classes <- ddata$class
outputs2 <- data.matrix(ddata[,c('output0','output1')])
softmax2 <- apply(outputs2, 1, function(x){exp(x[2])/sum(exp(x))})

##
plan(sequential)
plan(multisession, workers=6)
probs1 <- rowMeans(samplesF(Y=cbind(class=1), X=outputs2, parmList=parmlist, inorder=F))


## These functions make sure to chose equally in case of tie
maxdraw <- function(x){ (if(x[1]==x[2]){-1}else{which.max(x)-1}) }
##
buildcm <- function(trueclasses, probs, um=diag(2)){
    if(is.null(dim(probs))){probs <- rbind(1-probs, probs)}
    choices <- apply(um %*% probs, 2, maxdraw)
    ##
    cm <- matrix(c(
        sum(trueclasses==0 & choices==0),
        sum(trueclasses==0 & choices==1),
        sum(trueclasses==1 & choices==0),
        sum(trueclasses==1 & choices==1)
    ), 2, 2) +
        matrix(c(
            sum(trueclasses==0 & choices==-1),
            sum(trueclasses==0 & choices==-1),
            sum(trueclasses==1 & choices==-1),
            sum(trueclasses==1 & choices==-1)
        ), 2, 2)/2
    cm
}
##
comparescores <- function(trueclasses, um, outputs, probs, softmaxs){
    c('standard'=sum(buildcm(trueclasses, outputs) * um),
      'mixed'=sum(buildcm(trueclasses, softmaxs, um) * um),
      'transducer'=sum(buildcm(trueclasses, probs, um) * um)
      )
}

ulist <- list(c(1,0,0,1),
                    c(1,-10,0,10),
                    c(1,-100,0,100),
                    c(10,0,-10,1),
                    c(100,0,-100,1))
##
umlist <- lapply(ulist, function(x){
        x <- x-min(x)
        x <- x/max(x)
        matrix(x,2,2, byrow=T)
    }
)
##
ulist2 <- unlist(ulist)
dim(ulist2) <- c(4,length(ulist))
ulist2 <- t(ulist2)

results1 <- t(sapply(umlist, function(um){
    comparescores(trueclasses=classes, um=um, outputs=t(outputs2), probs=probs1, softmaxs=softmax2)/length(classes)}))
##
rresults <- round(results1,3)
cbind(ulist2, rresults,
      'rel_diff_std'=round(100*apply(rresults,1,function(x){diff(x[c(1,3)])/abs(x[1])}),1),
      'rel_diff_mix'=round(100*apply(rresults,1,function(x){diff(x[c(2,3)])/abs(x[2])}),1))


##                        standard mixed transducer rel_diff_std rel_diff_mix
## [1,]   1    0    0   1    0.959 0.959      0.961          0.2          0.2
## [2,]   1  -10    0  10    0.576 0.581      0.582          1.0          0.2
## [3,]   1 -100    0 100    0.536 0.548      0.548          2.2          0.0
## [4,]  10    0  -10   1    0.931 0.946      0.955          2.6          1.0
## [5,] 100    0 -100   1    0.928 0.950      0.955          2.9          0.5

#### with original utility values
##                        standard  mixed transducer rel_diff_std rel_diff_mix
## [1,]   1    0    0   1    0.959  0.959      0.961          0.2          0.2
## [2,]   1  -10    0  10    1.518  1.628      1.646          8.4          1.1
## [3,]   1 -100    0 100    7.237  9.526      9.573         32.3          0.5
## [4,]  10    0  -10   1    8.628  8.914      9.091          5.4          2.0
## [5,] 100    0 -100   1   85.584 89.941     90.914          6.2          1.1





## Calculation of utility yields on demonstration set
## for uniform distribution of utility matrices

id <- diag(2)
utp <- matrix(c(1,0,0,0),2,2)
utn <- matrix(c(0,0,0,1),2,2)
ufp <- matrix(c(0,0,1,0),2,2)
ufn <- matrix(c(0,1,0,0),2,2)
xy2um <- function(ab,ab2=NULL,norm=TRUE){
    if(length(ab)==1){ab <- c(ab,ab2)}
        um <- id + (ab[1]<0)*ab[1]*utn - (ab[1]>0)*ab[1]*utp +
            (ab[2]>0)*ab[2]*ufp - (ab[2]<0)*ab[2]*ufn
    if(norm){
        um <- um-min(um)
        um <- um/max(um)
    }
    um
}

nn <- 10^4
nn2 <- nn#3*10^3
##
lxy <- runif(2*round(nn*6/3),-1,1)
##
dim(lxy) <- c(round(nn*6/3),2)
lxy <- lxy[lxy[,2]<lxy[,1]+1 & lxy[,2]>lxy[,1]-1 & abs(lxy[,1])<=1 & abs(lxy[,2])<=1,][1:nn,]
##
lut <- apply(lxy,1,xy2um)
dim(lut) <- c(2,2,nn)

cmstandard <- buildcm(classes, t(outputs2))
##
allscores <- apply(lut, 3, function(um){
    c(sum(cmstandard * um),
      sum(buildcm(classes, softmax2, um) * um),
      sum(buildcm(classes, probs1, um) * um)
      )
})/length(classes)
## allscores <- (foreach(i=1:nn, .combine=cbind, .inorder=F)%dopar%{
##     um <- lut[,,i]
##     c(sum(cmstandard * um),
##       sum(buildcm(classes, probs1, um) * um),
##       sum(buildcm(classes, outputs1, um) * um)
##       )
## })/length(classes)
rownames(allscores) <- c('standard', 'mixed', 'transducer')

rowMeans(allscores)
#### skip=4
##  standard      mixed transducer 
## 0.7561422  0.7590596  0.7605033 
#### skip=8
##   standard      mixed transducer 
##  0.7545872  0.7573737  0.7587412 

pdff('../CNN_transducer_gainsx', asp=1)
## tplot(x=allscores[1,1:nn2], y=allscores[2,1:nn2]-allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5)
tplot(x=log10(allscores[1,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, asp=1,
      xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
      xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
      yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
      ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
      xlab='utility yield with standard method',
      ylab='utility yield with transducer & utility maximization')
abline(0,1, col=paste0(palette()[2],'88'), lwd=2, lty=1)
## tplot(x=log10(allscores[2,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, col=2,
##       xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       xlab='utility yield from mixed method',
##       ylab='utility yield from transducer & utility maximization')
## abline(0,1, col=paste0(palette()[4],'88'), lwd=2, lty=1)
dev.off()



## pdff('CNN_transducer_gains_relative')
## ## tplot(x=allscores[1,1:nn2], y=allscores[2,1:nn2]-allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5)
## tplot(x=allscores[1,1:nn2], y=100*(allscores[3,1:nn2]-allscores[1,1:nn2])/allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5,
##       xlab='utility yield from standard classification',
##       ylab='increase in utility yield/%')
## tplot(x=allscores[2,1:nn2], y=100*(allscores[3,1:nn2]-allscores[2,1:nn2])/allscores[2,1:nn2], type='p', pch=16, cex=1, alpha=0.5, col=2,
##       xlab='utility yield from mixed method',
##       ylab='increase in utility yield/%')
## dev.off()
## ##
## pdff('CNN_transducer_gains')
## ## tplot(x=allscores[1,1:nn2], y=allscores[2,1:nn2]-allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5)
## tplot(x=log10(allscores[1,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25,
##       xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       xlab='utility yield from standard classification',
##       ylab='utility yield from transducer & utility maximization')
## abline(0,1, col=paste0(palette()[4],'88'), lwd=2, lty=1)
## tplot(x=log10(allscores[2,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, col=2,
##       xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       xlab='utility yield from mixed method',
##       ylab='utility yield from transducer & utility maximization')
## abline(0,1, col=paste0(palette()[4],'88'), lwd=2, lty=1)
## dev.off()

#########################################################
## Calculation of utility yields on demonstration set
## altered base rates, generative mode
#########################################################
discardp <- sample(which(ddata$class==0), size=sum(ddata$class==0)-round(sum(ddata$class==1)/2))
classesb <- classes[-discardp]
outputs2b <- outputs2[-discardp,]
## softmax2 <- apply(outputs2, 1, function(x){exp(x[2])/sum(exp(x))})

##
probs1b <- probs1[-discardp]

probsinv1 <- rowMeans(samplesF(Y=outputs2b, X=cbind(class=1), parmList=parmlist, inorder=F))
##
probsinv0 <- rowMeans(samplesF(Y=outputs2b, X=cbind(class=0), parmList=parmlist, inorder=F))

baserates <- c('0'=sum(classesb==0), '1'=sum(classesb==1))/length(classesb)
##
bayesprobs1 <- probsinv1*baserates['1']/(probsinv0*baserates['0'] + probsinv1*baserates['1'])

## These functions make sure to chose equally in case of tie
maxdraw <- function(x){ (if(x[1]==x[2]){-1}else{which.max(x)-1}) }
##
buildcm <- function(trueclasses, probs, um=diag(2)){
    if(is.null(dim(probs))){probs <- rbind(1-probs, probs)}
    choices <- apply(um %*% probs, 2, maxdraw)
    ##
    cm <- matrix(c(
        sum(trueclasses==0 & choices==0),
        sum(trueclasses==0 & choices==1),
        sum(trueclasses==1 & choices==0),
        sum(trueclasses==1 & choices==1)
    ), 2, 2) +
        matrix(c(
            sum(trueclasses==0 & choices==-1),
            sum(trueclasses==0 & choices==-1),
            sum(trueclasses==1 & choices==-1),
            sum(trueclasses==1 & choices==-1)
        ), 2, 2)/2
    cm
}
##
comparescores <- function(trueclasses, um, outputs, probs, bayesprobs){
    c('standard'=sum(buildcm(trueclasses, outputs) * um),
      'discr'=sum(buildcm(trueclasses, probs, um) * um),
      'gener'=sum(buildcm(trueclasses, bayesprobs, um) * um)
      )
}

ulist <- list(c(1,0,0,1),
                    c(1,-10,0,10),
                    c(1,-100,0,100),
                    c(10,0,-10,1),
                    c(100,0,-100,1))
##
umlist <- lapply(ulist, function(x){
        x <- x-min(x)
        x <- x/max(x)
        matrix(x,2,2, byrow=T)
    }
)
##
ulist2 <- unlist(ulist)
dim(ulist2) <- c(4,length(ulist))
ulist2 <- t(ulist2)

results1 <- t(sapply(umlist, function(um){
    comparescores(trueclasses=classesb, um=um, outputs=t(outputs2b), probs=probs1b, bayesprobs=bayesprobs1)/length(classesb)}))
##
rresults <- round(results1,3)
cbind(ulist2, rresults,
      'rel_diff_std'=round(100*apply(rresults,1,function(x){diff(x[c(1,3)])/abs(x[1])}),1),
      'rel_diff_discr'=round(100*apply(rresults,1,function(x){diff(x[c(2,3)])/abs(x[2])}),1))

#### skip=8
## > +                        standard discr gener rel_diff_std rel_diff_discr
## [1,]   1    0    0   1    0.888 0.838 0.937          5.5           11.8
## [2,]   1  -10    0  10    0.749 0.824 0.838         11.9            1.7
## [3,]   1 -100    0 100    0.735 0.832 0.833         13.3            0.1
## [4,]  10    0  -10   1    0.683 0.667 0.680         -0.4            1.9
## [5,] 100    0 -100   1    0.657 0.667 0.667          1.5            0.0
#### skip=4
##                        standard discr gener rel_diff_std rel_diff_discr
## [1,]   1    0    0   1    0.896 0.849 0.941          5.0           10.8
## [2,]   1  -10    0  10    0.750 0.824 0.839         11.9            1.8
## [3,]   1 -100    0 100    0.735 0.832 0.833         13.3            0.1
## [4,]  10    0  -10   1    0.691 0.667 0.690         -0.1            3.4
## [5,] 100    0 -100   1    0.665 0.667 0.667          0.3            0.0


## Calculation of utility yields on demonstration set
## for uniform distribution of utility matrices

## id <- diag(2)
## utp <- matrix(c(1,0,0,0),2,2)
## utn <- matrix(c(0,0,0,1),2,2)
## ufp <- matrix(c(0,0,1,0),2,2)
## ufn <- matrix(c(0,1,0,0),2,2)
## xy2um <- function(ab,ab2=NULL,norm=TRUE){
##     if(length(ab)==1){ab <- c(ab,ab2)}
##         um <- id + (ab[1]<0)*ab[1]*utn - (ab[1]>0)*ab[1]*utp +
##             (ab[2]>0)*ab[2]*ufp - (ab[2]<0)*ab[2]*ufn
##     if(norm){
##         um <- um-min(um)
##         um <- um/max(um)
##     }
##     um
## }

## nn <- 10^4
## nn2 <- nn#3*10^3
## ##
## lxy <- runif(2*round(nn*6/3),-1,1)
## ##
## dim(lxy) <- c(round(nn*6/3),2)
## lxy <- lxy[lxy[,2]<lxy[,1]+1 & lxy[,2]>lxy[,1]-1 & abs(lxy[,1])<=1 & abs(lxy[,2])<=1,][1:nn,]
## ##
## lut <- apply(lxy,1,xy2um)
## dim(lut) <- c(2,2,nn)

cmstandard <- buildcm(classesb, t(outputs2b))
##
allscoresb <- apply(lut, 3, function(um){
    c(sum(cmstandard * um),
      sum(buildcm(classesb, probs1b, um) * um),
      sum(buildcm(classesb, bayesprobs1, um) * um)
      )
})/length(classesb)
## allscores <- (foreach(i=1:nn, .combine=cbind, .inorder=F)%dopar%{
##     um <- lut[,,i]
##     c(sum(cmstandard * um),
##       sum(buildcm(classes, probs1, um) * um),
##       sum(buildcm(classes, outputs1, um) * um)
##       )
## })/length(classes)
rownames(allscoresb) <- c('standard', 'discr', 'gener')

rowMeans(allscoresb)
#### skip=8
##  standard     discr     gener 
## 0.7141566 0.7079744 0.7451940 

#### skip=4
##  standard     discr     gener 
## 0.7190232 0.7112459 0.7479897 

pdff('../CNN_transducer_gains_generx', asp=1)
## tplot(x=allscores[1,1:nn2], y=allscores[2,1:nn2]-allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5)
tplot(x=log10(allscoresb[1,1:nn2]), y=log10(allscoresb[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, asp=1,
      xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
      xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
      yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
      ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
      xlab='utility yield with standard method',
      ylab='utility yield with transducer, generative mode')
abline(0,1, col=paste0(palette()[2],'88'), lwd=2, lty=1)
## tplot(x=log10(allscores[2,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, col=2,
##       xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       xlab='utility yield from mixed method',
##       ylab='utility yield from transducer & utility maximization')
## abline(0,1, col=paste0(palette()[4],'88'), lwd=2, lty=1)
dev.off()
##
pdff('../CNN_transducer_gains_gener_vs_discrx', asp=1)
## tplot(x=allscores[1,1:nn2], y=allscores[2,1:nn2]-allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5)
tplot(x=log10(allscoresb[2,1:nn2]), y=log10(allscoresb[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, asp=1, col=3,
      xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
      xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
      yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
      ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
      xlab='utility yield with transducer, discriminative mode',
      ylab='utility yield with transducer, generative mode')
abline(0,1, col=paste0(palette()[2],'88'), lwd=2, lty=1)
## tplot(x=log10(allscores[2,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, col=2,
##       xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       xlab='utility yield from mixed method',
##       ylab='utility yield from transducer & utility maximization')
## abline(0,1, col=paste0(palette()[4],'88'), lwd=2, lty=1)
dev.off()


## pdff('CNN_transducergen_gains_relative')
## ## tplot(x=allscores[1,1:nn2], y=allscores[2,1:nn2]-allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5)
## tplot(x=allscores[1,1:nn2], y=100*(allscores[3,1:nn2]-allscores[1,1:nn2])/allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5,
##       xlab='utility yield from standard classification',
##       ylab='increase in utility yield/%')
## tplot(x=allscores[2,1:nn2], y=100*(allscores[3,1:nn2]-allscores[2,1:nn2])/allscores[2,1:nn2], type='p', pch=16, cex=1, alpha=0.5, col=2,
##       xlab='utility yield from mixed method',
##       ylab='increase in utility yield/%')
## dev.off()

## pdff('CNN_transducergen_gains', asp=1)
## ## tplot(x=allscores[1,1:nn2], y=allscores[2,1:nn2]-allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5)
## tplot(x=log10(allscores[1,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25,
##       xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       xlab='utility yield with standard method',
##       ylab='utility yield with transducer & utility maximization')
## abline(0,1, col=paste0(palette()[2],'88'), lwd=2, lty=1)
## ## tplot(x=log10(allscores[2,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, col=2,
## ##       xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
## ##       xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
## ##       yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
## ##       ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
## ##       xlab='utility yield from mixed method',
## ##       ylab='utility yield from transducer & utility maximization')
## ## abline(0,1, col=paste0(palette()[4],'88'), lwd=2, lty=1)
## dev.off()


#########################################################
## 
#########################################################












xgrid <- seq(0, 1, length.out=256)
ygrid <- X2Y[[outputcov]](xgrid)
##
vpoints <- cbind(ygrid)
colnames(vpoints) <- outputcov
##
opgrid <- samplesF(Y=cbind(class=1), X=vpoints, parmList=parmlist, inorder=F)
##
qgrid <- apply(opgrid,1,function(x){quantile(x, c(1,7)/8)})
##

tplot(x=xgrid, y=cbind(rowMeans(opgrid), 1- rowMeans(opgrid)), xlab='output',
##      ylab=expression(p~group('(',class~output,')')),
      ylab=bquote('P'~group('(','class', '.')~group('|', ' output',')')),
      mar=c(4.5,5.5,1,1),
      ylim=c(0,1), lwd=3, family='Palatino', asp=1)
legend(x=0.25,y=1.05, c('class 1', 'class 0'), lty=c(1,2), col=c(1,2), lwd=3, bty='n', cex=1.25)
##
polygon(x=c(xgrid,rev(xgrid)), y=c(qgrid[1,],rev(qgrid[2,])), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=1-c(qgrid[1,],rev(qgrid[2,])), col=paste0(palette()[2],'40'), border=NA)



## legend('topleft', legend=c(
##                        paste0(paste0(rownames(qgrid),collapse='\u2013'), ' uncertainty')
##                    ),
##        lty=c('solid'), lwd=c(10),
##        col=paste0(palette()[1],c('40')),
##        bty='n', cex=1.25)


#########################################################
## Check of Kjetil's calculations test set 2
#########################################################
kresults <- fread('CNN_bayesian_prob_test2.csv', sep=',')
colnames(kresults)[2:3] <- realCovs

ptest1 <- samplesF(Y=cbind(class=0), X=data.matrix(kresults[,..realCovs]), parmList=parmlist, inorder=F)

mptest1 <- rowMeans(ptest1)

tplot(x=kresults$prob_0, y=mptest1, type='p', xlab='Kjetil',ylab='Luca', pch='+')

discre <- (abs(mptest1-kresults[,prob_0])/mptest1)
which.max(discre)

summary(discre)*100

#########################################################
## Comparison of utility scores with & without Bayesian augmentation
#########################################################
classseq <- kresults$class

kresults$sigmoid <- apply(kresults[,..realCovs],1,function(x){exp(x[1])/sum(exp(x))})

kscores <- fread('scores_utility_matrix.csv', sep=',')

maxdraw <- function(x){ (if(x[1]==x[2]){-1}else{which.max(x)-1}) }
##
uscore <- function(trueclass, prob, um=diag(2), umchoice=um){
    choices <- apply(umchoice %*% rbind(prob,1-prob), 2, maxdraw)
    ##
    cm <- matrix(c(
        sum(trueclass==0 & choices==0),
        sum(trueclass==0 & choices==1),
        sum(trueclass==1 & choices==0),
        sum(trueclass==1 & choices==1)
    ), 2, 2) +
        matrix(c(
        sum(trueclass==0 & choices==-1),
        sum(trueclass==0 & choices==-1),
        sum(trueclass==1 & choices==-1),
        sum(trueclass==1 & choices==-1)
    ), 2, 2)/2
    sum(um * cm)
}
##
comparescores <- function(trueclass, um, output, ...){
    c(
        uscore(trueclass, output, um=um, umchoice=diag(2)),
        sapply(list(...),function(x){uscore(trueclass, x, um)})
    )
}

## umlist <- lapply(
##     list(c(1,0,0,1),
##          c(4,0,0,1),
##          c(32,0,0,1),
##          c(1,0,0,4),
##          c(1,0,0,32),
##          c(1,-1,0,1),
##          c(1,-4,0,1),
##          c(1,-32,0,1),
##          c(1,0,-1,1),
##          c(1,0,-4,1),
##          c(1,0,-32,1)),
##     function(x){
##         x <- x-min(x)
##         x <- x/max(x)
##         matrix(x,2,2)
##     }
## )

umlist <- lapply(
    list(c(1,0,0,1),
         c(1,0,-1,1),
         c(1,-1,0,1),
         c(10,0,0,1),
         c(1,0,0,10),
         c(1,0,-10,1),
         c(1,-10,0,1),
         c(5,0,0,1),
         c(100,0,0,1),
         c(1,0,0,5),
         c(1,0,0,100),
         c(1,-5,0,1),
         c(1,0,-5,1),
         c(1,-100,0,1),
         c(1,0,-100,1)),
    function(x){
        x <- x-min(x)
        x <- x/max(x)
        matrix(x,2,2)
    }
)

## results1b <- t(sapply(umlist, function(um){
##     comparescores2(kresults$class, um=um, kresults$sigmoid, kresults$sigmoid, mptest1)
## }))/length(kresults$class)

## using the full MCMC results
results1 <- t(sapply(umlist, function(um){
    comparescores(kresults$class, kresults$sigmoid, mptest1, um=um)
}))/length(kresults$class)
colnames(results1) <- c('raw', 'out', 'bayes')
##
cbind(results1,t(apply(results1,1,function(x){
    mx <- max(x)
    c(x-mx)/mx*100
}
)))

## using a reduced set of coefficients
results2 <- t(sapply(umlist, function(um){
    comparescores(kresults$class, kresults$sigmoid, kresults$prob_0, um=um)
}))/length(kresults$class)
colnames(results2) <- c('raw', 'out', 'bayes')
##
cbind(results2,t(apply(results2,1,function(x){
    mx <- max(x)
    c(x-mx)/mx*100
}
)))

cbind(results[,c(1,3)],t(apply(results[,c(1,3)],1,function(x){
    mx <- max(x)
    c(x-mx)/mx*100
}
)))

cbind(kscores[,c('score_CNN_output','score_CNN_bayesian')],t(apply(kscores[,c('score_RF_output','score_RF_bayesian')],1,function(x){
    mx <- max(x)
    c(x-mx)/mx*100
}
)))

cbind(results[,'raw'], kscores[,score_RF_output], results[,'raw']- kscores[,score_RF_output])*length(trueclass)

cbind(results[,'bayes'], kscores[,score_RF_bayesian], results[,'bayes']- kscores[,score_RF_bayesian])*length(trueclass)



#########################################################
## Comparison of utility scores with & without Bayesian augmentation
## for a large number of possible utility matrices
#########################################################

lo <- 128
ss <- seq(-1,1,length.out=lo)
xy <- cbind(rep(ss,lo), rep(ss,each=lo))
#xy <- xy[xy[,2]<=xy[,1]+1 & xy[,2]>=xy[,1]-1,]
##tplot(x=xy[,1], y=xy[,2], type='p')
##
id <- diag(2)
utp <- matrix(c(1,0,0,0),2,2)
utn <- matrix(c(0,0,0,1),2,2)
ufp <- matrix(c(1,0,1,0),2,2)
ufn <- matrix(c(1,1,0,0),2,2)
##
resultsl <- t(future_apply(xy, 1, function(coo){
    if(coo[2]<=coo[1]+1 & coo[2]>=coo[1]-1){
    um <- id + (coo[1]<0)*coo[1]*utn - (coo[1]>0)*coo[1]*utp +
        (coo[2]>0)*coo[2]*ufp - (coo[2]<0)*coo[2]*ufn
    ##
    comparescores(classseq, um=um, kresults$sigmoid, kresults$sigmoid, mptest1)
    }else{rep(NA,3)}
}))
colnames(resultsl) <- c('standard', 'output_as_prob', 'bayes')

reldiff <- round(100*(resultsl[,3] - resultsl[,1])/resultsl[,1],1)
summary(reldiff)
maxr <- max(reldiff,na.rm=T)
lreldiff <- plogis(reldiff)
##
rg <- max(abs(lreldiff), na.rm=T)#quantile(lreldiff,7/8, na.rm=T)
colr <- colorRamp(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3', '#2166ac'), space='rgb', interpolate='linear')
cex <- 0.65
nlevels <- 16
levels <- seq(0.5,plogis(maxr),length.out=nlevels)
baseseq <- qlogis(levels)
pdff('CNN_improv_umspace')
dim(lreldiff) <- c(lo,lo)
filled.contour(x=ss, y=ss, z=lreldiff, zlim=c(0.5,1), color.palette=function(n){rgb(colr(seq(0.5,1,length.out=n))/255)}, levels=levels, asp=1,frame=F,axes=F
               ,key.axes = axis(4, at=levels,
                               labels=round(baseseq,1),
                                asp=1)
               )
dev.off()

reldiff <- round(100*(resultsl[,3] - resultsl[,1])/resultsl[,1],1)
summary(reldiff)
rg <- range(c(0,reldiff))
groupp <- reldiff>= 0
groupn <- reldiff < 0
##
pdff('improvement_CNN')
tplot(x=resultsl[groupp,1]/length(classseq), y=reldiff[groupp], type='p', pch=16,cex=1, col=1,
      xlab='standard method', ylab='% improvement with bayes augmentation', xlim=0:1, ylim=rg)
if(sum(groupn)>0){
    tplot(x=resultsl[groupn,1]/length(classseq), y=reldiff[groupn], type='p', pch=16,cex=1, col=2, add=T)
}
dev.off()
sum(groupp)/length(reldiff)




## nn <- 10^4
## basetp <- matrix(c(1,0,0,0),2,2)
## basetn <- matrix(c(0,0,0,1),2,2)
## basefpfn <- array(c(1/2,1/2,0,0,  0,0,1/2,1/2),dim=c(2,2,2))
## sides <- sample(1:2,nn,replace=T)
## convpoints <- LaplacesDemon::rdirichlet(n=nn, alpha=rep(1,3))
## lum <- future_lapply(1:nn,function(i){
##     temp <- basetp * convpoints[i,1] + basetn * convpoints[i,2] +
##         basefpfn[,,sides[i]] * convpoints[i,3]
##     ## temp <- temp - min(temp)
##     temp <- temp/max(temp)
## })
## ##dim(lum) <- c(2,2,nn)

## resultsl <- t(future_sapply(lum, function(um){
##     comparescores(kresults$class, kresults$sigmoid, mptest1, um=um)
## }))/length(kresults$class)
## colnames(resultsl) <- c('raw', 'out', 'bayes')

## reldiff <- (resultsl[,3]-resultsl[,1])/resultsl[,1]*100
## summary(reldiff)
## rg <- range(c(0,reldiff))
## groupp <- reldiff>= 0
## groupn <- reldiff < 0
## ##
## pdff('improvement_CNN')
## tplot(x=resultsl[groupp,1], y=reldiff[groupp], type='p', pch='.', col=1,
##       xlab='standard method', ylab='% improvement with bayes augmentation', xlim=0:1, ylim=rg)
## if(sum(groupn)>0){
##     tplot(x=resultsl[groupn,1], y=reldiff[groupn], type='p', pch='.', col=2, add=T)
## }
## dev.off()


#########################################################
## Comparison of scores for a test set with different
## proportions of classes
#########################################################
classseq <- kresults$class

csubset <- sort(c(which(classseq==1),
                  sample(which(classseq==0), sum(classseq==1))))

probinv0 <- samplesF(X=cbind(class=0), Y=data.matrix(kresults[csubset,..realCovs]), parmList=parmlist, inorder=F)
probinv1 <- samplesF(X=cbind(class=1), Y=data.matrix(kresults[csubset,..realCovs]), parmList=parmlist, inorder=F)

mprobinv0 <- rowMeans(probinv0)
mprobinv1 <- rowMeans(probinv1)

probinvbase <- mprobinv0/(mprobinv0+mprobinv1)

## using the full MCMC results
resultsbase <- t(sapply(umlist, function(um){
    comparescores(classseq[csubset], um=um,
                  kresults$sigmoid[csubset],
                  kresults$sigmoid[csubset],
                  mptest1[csubset],
                  probinvbase)
}))
colnames(resultsbase) <- c('standard', 'output_as_prob', 'bayes', 'inverse_bayes')

##
cbind(resultsbase[,],t(apply(resultsbase[,],1,function(x){
    mx <- max(x)
    c(x-mx)/mx*100
}
)))

## using the full MCMC results
results1 <- t(sapply(umlist, function(um){
    comparescores(kresults$class[csubset], kresults$sigmoid[csubset], probinvbase, um=um)
}))
colnames(results1) <- c('raw', 'out', 'bayes')
##
cbind(results1[,c(1,3)],t(apply(results1[,c(1,3)],1,function(x){
    mx <- max(x)
    c(x-mx)/mx*100
}
)))

## using a reduced set of coefficients
results2 <- t(sapply(umlist, function(um){
    comparescores(kresults$class, kresults$sigmoid, kresults$prob_0, um=um)
}))/length(kresults$class)
colnames(results2) <- c('raw', 'out', 'bayes')
##
cbind(results2,t(apply(results2,1,function(x){
    mx <- max(x)
    c(x-mx)/mx*100
}
)))




#########################################################
## Check of Kjetil's calculations test set 1
#########################################################
kresults <- fread('CNN_direct_prob.csv', sep=',')

ptest1 <- samplesF(Y=cbind(class=0), X=data.matrix(kresults[,..realCovs]), parmList=parmlist, inorder=F)

mptest1 <- rowMeans(ptest1)

discre <- (abs(mptest1-kresults[,direct_prob_0])/mptest1)
which.max(discre)

summary(discre)*100
#########################################################



xr <- range(unlist(Xrange))
cseq <- seq(xr[1], xr[2], length.out=128)
##
vpoints <- cbind(prediction0=rep(cseq, length(cseq)), prediction1=rep(cseq, each=length(cseq)))

system.time(opgrid <- samplesF(Y=cbind(class=1), X=vpoints, parmList=parmlist, inorder=F))
##

## system.time(opgrid2 <- samplesFp(Y=cbind(class=1), X=vpoints, batchsize=1024, parmList=parmlist, inorder=F))


mpgrid <- rowMeans(opgrid)
dim(mpgrid) <- rep(length(cseq), 2)

fig <- plot_ly(z=t(mpgrid), x=cseq, y=cseq, cmin=0, cmax=1)
fig <- fig %>% add_surface(colors='Reds')
fig <- fig %>% layout(scene = list(xaxis = list(title = "output 0"), yaxis = list(title = "output 1"), zaxis = list(title = "probability of class 1", range=c(0,1)), camera = list(projection = list(type = 'orthographic'))), title='original')
fig
orca(fig, 'CNNprobability_vs_output1.pdf')


####
iqrgrid <- apply(pgrid,1,IQR)
dim(iqrgrid) <- rep(length(cseq), 2)

fig <- plot_ly(z=t(iqrgrid), x=cseq, y=cseq, cmin=0, cmax=1)
fig <- fig %>% add_surface(colors='Reds')
fig <- fig %>% layout(scene = list(xaxis = list(title = "output 0"), yaxis = list(title = "output 1"), zaxis = list(title = "interquartile range around probability of class 1", range=c(0,1))))
fig

####
q1grid <- apply(pgrid,1,function(x){quantile(x, 1/8)})
dim(q1grid) <- rep(length(cseq), 2)
q2grid <- apply(pgrid,1,function(x){quantile(x, 7/8)})
dim(q2grid) <- rep(length(cseq), 2)

fig <- plot_ly(x=cseq, y=cseq, cmin=0, cmax=1)
fig <- fig %>% add_surface(z = t(mpgrid), opacity = 0.5, colorscale='Blues')
fig <- fig %>% add_surface(z = t(q1grid), opacity = 0.5, colorscale='Reds')
fig <- fig %>% add_surface(z = t(q2grid), opacity = 0.5, colorscale='Reds')
fig <- fig %>% layout(scene = list(xaxis = list(title = "output 0"), yaxis = list(title = "output 1"), zaxis = list(title = "probability of class 1", range=c(0,1))))
fig


## Put all parameters into one list
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
##
qorder <- order(c(oneparmlist$q), decreasing=FALSE)
shortparmlist <- list(
    q=oneparmlist$q[,qorder, drop=F]/sum(oneparmlist$q[,qorder]),
    meanR=oneparmlist$meanR[,,qorder, drop=F],
    tauR=oneparmlist$tauR[,,qorder, drop=F],
    probI=oneparmlist$probI[,,qorder, drop=F],
    sizeI=oneparmlist$sizeI[,,qorder, drop=F],
    probB=oneparmlist$probB[,,qorder, drop=F]
)
fwrite(data.table(w=c(shortparmlist$q),
               p=c(shortparmlist$probB),
               mu0=c(shortparmlist$meanR[,'prediction0',]),
               sigma0=1/sqrt(c(shortparmlist$tauR[,'prediction0',])),
               mu1=c(shortparmlist$meanR[,'prediction1',]),
               sigma1=1/sqrt(c(shortparmlist$tauR[,'prediction1',]))
               ), 'CNN_probabilityfunction_full.csv', sep=',')


qorder <- order(c(oneparmlist$q), decreasing=TRUE)
#qcumsum <- cumsum(qlist[qorder])
## qorder <- order(c(oneparmList$q), decreasing=TRUE)
## ##
## maxclusters <- 2^10
## ##
## fq <- parmList$q[qorder[1:maxclusters]]
## fq <- fq/sum(fq)
## fmeanR <- oneparmlist$meanR[1,,qorder[1:maxclusters]]
## ftauR <- oneparmlist$tauR[1,,qorder[1:maxclusters]]
## fprobB <- oneparmlist$probB[1,,qorder[1:maxclusters]]
##
## qselect <- c(oneparmlist$q)[qorder] > 2^-64
qselect <- 1:(2^17)#(which(cumsum(c(oneparmlist$q)[qorder]) < 1-2^-16)[1])
maxclusters <- round(length(oneparmlist$q))
cincluded <- qorder[qselect]
##
shortparmlist <- list(
    q=oneparmlist$q[,cincluded, drop=F]/sum(oneparmlist$q[,cincluded]),
    meanR=oneparmlist$meanR[,,cincluded, drop=F],
    tauR=oneparmlist$tauR[,,cincluded, drop=F],
    probI=oneparmlist$probI[,,cincluded, drop=F],
    sizeI=oneparmlist$sizeI[,,cincluded, drop=F],
    probB=oneparmlist$probB[,,cincluded, drop=F]
)

##
## shortpgrid <- samplesF(Y=cbind(class=1), X=vpoints, parmList=shortparmlist, inorder=F)

shortpgrid <- foreach(i=1:(nrow(vpoints)/128), .combine=c)%dopar%{
    samplesF(Y=cbind(class=1), X=vpoints[(i-1)*128+(1:128),], parmList=shortparmlist, inorder=F)
}
dim(shortpgrid) <- rep(length(cseq), 2)


##
##
fig <- plot_ly(z=t(shortpgrid), x=cseq, y=cseq, cmin=0, cmax=1)
fig <- fig %>% add_surface(colors='Blues')
fig <- fig %>% layout(scene = list(xaxis = list(title = "output 0"), yaxis = list(title = "output 1"), zaxis = list(title = "probability of class 1", range=c(0,1))), title='approx')
fig




fwrite(data.table(w=c(shortparmlist$q),
               p=c(shortparmlist$probB[,1,]),
               mu0=c(shortparmlist$meanR[,1,]),
               sigma0=1/sqrt(c(shortparmlist$tauR[,1,])),
               mu1=c(shortparmlist$meanR[,2,]),
               sigma1=1/sqrt(c(shortparmlist$tauR[,2,]))
               ), 'CNN_probabilityfunction2.csv', sep=',')



##
xgrid <- seq(0, 1, length.out=256)
ygrid <- X2Y[['prediction_lnodds']](xgrid)
##
vpoints <- cbind(prediction_lnodds=ygrid)
##
pgrid <- samplesF(Y=cbind(class=1), X=vpoints, parmList=shortparmlist, inorder=F)
##
## qgrid <- apply(pgrid,1,function(x){quantile(x, c(1,7)/8)})
##
tplot(x=xgrid, y=1-pgrid, xlab='RF % output', ylab='probability of class 1', ylim=c(0,1))
## polygon(x=c(xgrid,rev(xgrid)), y=c(qgrid[1,],rev(qgrid[2,])), col=paste0(palette()[1],'40'), border=NA)
## legend('topleft', legend=c(
##                        paste0(paste0(rownames(qgrid),collapse='\u2013'), ' uncertainty')
##                    ),
##        lty=c('solid'), lwd=c(10),
##        col=paste0(palette()[1],c('40')),
##        bty='n', cex=1.25)



pgrid0 <- samplesF(X=cbind(class=0), Y=vpoints, parmList=parmlist, inorder=F)*Xjacobian[['prediction_lnodds']](xgrid)
##
qgrid0 <- apply(pgrid0,1,function(x){quantile(x, c(1,7)/8)})
##
pgrid1 <- samplesF(X=cbind(class=1), Y=vpoints, parmList=parmlist, inorder=F)*Xjacobian[['prediction_lnodds']](xgrid)
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


fwrite(data.table(w=c(shortparmlist$q),
               p=c(shortparmlist$probB),
               mu=c(shortparmlist$meanR),
               sigma=1/sqrt(c(shortparmlist$tauR))
               ), 'testRFfunction.csv', sep=',')





## Allocation period 2022.1 (2022-04-01 - 2022-10-01)
## Accounting updated 2022-04-29 17:33:20
## ============================================
## Account                            Cpu hours
## --------------------------------------------
## nn8050k  Quota (pri)                14000.00
## nn8050k  Quota (nonpri)                 0.00
## nn8050k  Used                        3257.52
## nn8050k  Running                        0.00
## nn8050k  Pending                       40.00
## nn8050k  Available                  10702.48
## ============================================











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

testj <- rowMeans(samplesF(Y=as.matrix(alldata[]), parmList=parmList, inorder=F))
testc <- rowMeans(samplesF(Y=as.matrix(alldata[,'class']), parmList=parmList, inorder=F))
testp <- rowMeans(samplesF(Y=as.matrix(alldata[,'prediction_lnodds']), parmList=parmList, inorder=F))
##
testcp <- rowMeans(samplesF(Y=as.matrix(alldata[,'class']), X=as.matrix(alldata[,'prediction_lnodds']), parmList=parmList, inorder=F))
testpc <- rowMeans(samplesF(X=as.matrix(alldata[,'class']), Y=as.matrix(alldata[,'prediction_lnodds']), parmList=parmList, inorder=F))

max(abs((testcp-testj/testp)/(testcp+testj/testp)))
max(abs((testpc-testj/testc)/(testpc+testj/testc)))




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

