## Author: PGL  Porta Mana
## Created: 2022-03-17T14:21:57+0100
## Last-Updated: 2022-05-27T23:54:22+0200
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
signif2 <- function(x, s){signif(x*2L, s)/2L}
round2 <- function(x, s){round(x*2L, s)/2L}
#### End custom setup ####

#########################################################
## Setup and read calibration data
#########################################################
set.seed(707)
baseversion <- '_rfcont_new'
maincov <- 'class'
outputcov <- 'transf_output1'
family <- 'Palatino'
saveinfofile <- 'rfcont_variateinfo.csv'
calibfile <- 'tmodCHEMBL205_predictions_RF_test1_calibration.csv'
demofile <- 'tmodCHEMBL205_predictions_RF_test2_demonstration.csv'
##
X2Y <- list(
    function(x){
        epsi <- 1 - 2^-10
        x <- 0.5 + (x-0.5)*epsi
        log(x/(1-x))
    }
)
Xjacobian <- list(
    function(x){
        epsi <- 1 - 2^-10
        x <- 0.5 + (x-0.5)*epsi
        epsi/(x*(1-x))
    }
)
Xrange <- list(c(0,1))
names(X2Y) <- names(Xjacobian) <- names(Xrange) <- outputcov

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
if(!exists('X2Y')){X2Y <- list()}
if(!exists('Xjacobian')){Xjacobian <- list()}
if(!exists('Xrange')){Xrange <- list()}
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
dirname <- '_rfcont_new-V2-D3589-K64-I1024'

npar <- 16
ntotal <- 1024*4
nskip <- 4
parmlist <- mcsamples2parmlist(
    foreach(i=1:npar, .combine=rbind)%dopar%{
        temp <- readRDS(paste0(dirname, '/_mcsamples-R_rfcont_new_',i,'_1-V2-D3589-K64-I1024.rds'))
        if(any(is.na(nrow(temp)+1-rev(seq(1,nrow(temp),by=nskip)[1:(ntotal/npar)])))){print('WARNING! not enough points')}
        temp[nrow(temp)+1-rev(seq(1,nrow(temp),by=nskip)[1:(ntotal/npar)]),]
}
)

#########################################################
## transducer curve p(c | y)
#########################################################

xgrid <- seq(0, 1, length.out=256)
ygrid <- X2Y[[outputcov]](xgrid)
##
vpoints <- cbind(ygrid)
colnames(vpoints) <- outputcov
##
plan(sequential)
plan(multisession, workers=6)
opgrid <- samplesF(Y=cbind(class=1), X=vpoints, parmList=parmlist, inorder=F)
##
qgrid <- apply(opgrid,1,function(x){quantile(x, c(1,7)/8)})
##

pdff('../transducer_curve_RFx', asp=1)
tplot(x=xgrid, y=cbind(rowMeans(opgrid), 1- rowMeans(opgrid)), xlab='output',
##      ylab=expression(p~group('(',class~output,')')),
      ylab=bquote('P'~group('(','class', '.')~group('|', ' output',')')),
      mar=c(4.5,5.5,1,1),
      ylim=c(0,1), lwd=3, family='Palatino', asp=1)
legend(x=0.25,y=1.05, c('class 1', 'class 0'), lty=c(1,2), col=c(1,2), lwd=3, bty='n', cex=1.5)
##
polygon(x=c(xgrid,rev(xgrid)), y=c(qgrid[1,],rev(qgrid[2,])), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=1-c(qgrid[1,],rev(qgrid[2,])), col=paste0(palette()[2],'40'), border=NA)
dev.off()

#########################################################
## transducer curve p(y | c)
#########################################################

xgrid <- seq(0, 1, length.out=256)
ygrid <- X2Y[[outputcov]](xgrid)
##
vpoints <- cbind(ygrid)
colnames(vpoints) <- outputcov

plan(sequential)
plan(multisession, workers=6)
py0grid <- samplesF(Y=vpoints, X=cbind(class=0), parmList=parmlist, inorder=F)
py1grid <- samplesF(Y=vpoints, X=cbind(class=1), parmList=parmlist, inorder=F)
##
q0grid <- apply(py0grid,1,function(x){quantile(x, c(1,7)/8)})
q1grid <- apply(py1grid,1,function(x){quantile(x, c(1,7)/8)})
##

pdff('../transducer_curve_RFx_inverse')
tplot(x=xgrid, y=cbind(rowMeans(py1grid), rowMeans(py0grid))*Xjacobian[[outputcov]](xgrid), xlab='output',
##      ylab=expression(p~group('(',class~output,')')),
      ylab=bquote('p'~group('(','output', '.')~group('|', ' class',')')),
      mar=c(4.5,5.5,1,1),
      ylim=c(0,3), lwd=3, family='Palatino')
legend('top', c('class 1', 'class 0'), lty=c(1,2), col=c(1,2), lwd=3, bty='n', cex=1.5)
##
polygon(x=c(xgrid,rev(xgrid)), y=c(q1grid[1,]*Xjacobian[[outputcov]](xgrid),rev(q1grid[2,]*Xjacobian[[outputcov]](xgrid))), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(q0grid[1,]*Xjacobian[[outputcov]](xgrid),rev(q0grid[2,]*Xjacobian[[outputcov]](xgrid))), col=paste0(palette()[2],'40'), border=NA)
dev.off()



## legend('topleft', legend=c(
##                        paste0(paste0(rownames(qgrid),collapse='\u2013'), ' uncertainty')
##                    ),
##        lty=c('solid'), lwd=c(10),
##        col=paste0(palette()[1],c('40')),
##        bty='n', cex=1.25)


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
               mu=c(shortparmlist$meanR),
               sigma=1/sqrt(c(shortparmlist$tauR))
               ), '_RF_transducer_parameters.csv', sep=',')


#########################################################
## Calculation of utility yields on demonstration set
#########################################################
ddata <- fread(demofile, sep=',')
classes <- ddata$class
outputs1 <- ddata$output1
transfoutputs1 <- cbind(X2Y[[1]](outputs1))
colnames(transfoutputs1) <- outputcov
##
plan(sequential)
plan(multisession, workers=6)
probs1 <- rowMeans(samplesF(Y=cbind(class=1), X=transfoutputs1, parmList=parmlist, inorder=F))


probsinv1 <- rowMeans(samplesF(Y=transfoutputs1, X=cbind(class=1), parmList=parmlist, inorder=F))*Xjacobian[[1]](outputs1)
##
probsinv0 <- rowMeans(samplesF(Y=transfoutputs1, X=cbind(class=0), parmList=parmlist, inorder=F))*Xjacobian[[1]](outputs1)

fwrite(cbind(ddata, 'RF_prob1'=probs1,
             'RF_invprob0'=probsinv0, 'RF_invprob1'=probsinv1),
       paste0('RFplus_',demofile), sep=',')


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
comparescores <- function(trueclasses, um, outputs, probs){
    c('standard'=sum(buildcm(trueclasses, outputs) * um),
      'mixed'=sum(buildcm(trueclasses, outputs, um) * um),
      'transducer'=sum(buildcm(trueclasses, probs, um) * um)
      )
}

ulist <- list(c(1,0,0,1),
              c(1,-10,0,10),
              c(1,-100,0,100),
              c(10,0,-10,1),
              c(100,0,-100,1),
              ##
              c(1,-10,-1,10),
              c(1,-100,-1,100),
              c(10,-1,-10,1),
              c(100,-1,-100,1)
              )
##
umlist <- c(lapply(ulist, function(x){ matrix(x,2,2, byrow=T) } ),
    lapply(ulist, function(x){
        x <- x-min(x)
        x <- x/max(x)
        matrix(x,2,2, byrow=T)
    }
) )
##
ulist2 <- unlist(umlist)
dim(ulist2) <- c(4,2*length(ulist))
ulist2 <- t(ulist2)

results1 <- t(sapply(umlist, function(um){
    comparescores(trueclasses=classes, um=um, outputs=outputs1, probs=probs1)/length(classes)}))
##
rresults <- round(results1,3)
##
options(width=160)
cbind(ulist2, rresults,
      'd_std'=round(apply(rresults,1,function(x){diff(x[c(1,3)])}),4),
      'd_mix'=round(apply(rresults,1,function(x){diff(x[c(2,3)])}),4),
      'rd%_std'=round(100*apply(rresults,1,function(x){diff(x[c(1,3)])/abs(x[1])}),2),
      'rd%_mix'=round(100*apply(rresults,1,function(x){diff(x[c(2,3)])/abs(x[2])}),2))

##                                         standard  mixed transducer d_std  d_mix rd%_std rd%_mix
##  [1,]   1.000    0.000    0.000   1.000    0.968  0.968      0.974 0.006  0.006    0.62    0.62
##  [2,]   1.000    0.000  -10.000  10.000    1.364  1.555      1.719 0.355  0.164   26.03   10.55
##  [3,]   1.000    0.000 -100.000 100.000    5.553  9.132      9.617 4.064  0.485   73.19    5.31
##  [4,]  10.000  -10.000    0.000   1.000    8.954  9.083      9.091 0.137  0.008    1.53    0.09
##  [5,] 100.000 -100.000    0.000   1.000   88.920 90.915     90.914 1.994 -0.001    2.24    0.00
##
##  [6,]   1.000   -1.000  -10.000  10.000    1.354  1.603      1.681 0.327  0.078   24.15    4.87
##  [7,]   1.000   -1.000 -100.000 100.000    5.543  8.433      9.509 3.966  1.076   71.55   12.76
##  [8,]  10.000  -10.000   -1.000   1.000    8.932  8.996      9.001 0.069  0.005    0.77    0.06
##  [9,] 100.000 -100.000   -1.000   1.000   88.898 90.824     90.823 1.925 -0.001    2.17    0.00
####
## [10,]   1.000    0.000    0.000   1.000    0.968  0.968      0.974 0.006  0.006    0.62    0.62
## [11,]   0.550    0.500    0.000   1.000    0.568  0.578      0.586 0.018  0.008    3.17    1.38
## [12,]   0.505    0.500    0.000   1.000    0.528  0.546      0.548 0.020  0.002    3.79    0.37
## [13,]   1.000    0.000    0.500   0.550    0.948  0.954      0.955 0.007  0.001    0.74    0.10
## [14,]   1.000    0.000    0.500   0.505    0.945  0.955      0.955 0.010  0.000    1.06    0.00
##
## [15,]   0.550    0.450    0.000   1.000    0.568  0.580      0.584 0.016  0.004    2.82    0.69
## [16,]   0.505    0.495    0.000   1.000    0.528  0.542      0.548 0.020  0.006    3.79    1.11
## [17,]   1.000    0.000    0.450   0.550    0.947  0.950      0.950 0.003  0.000    0.32    0.00
## [18,]   1.000    0.000    0.495   0.505    0.944  0.954      0.954 0.010  0.000    1.06    0.00


t(cbind(ulist2, signif(rresults[,c(1,3)],3), round(100*apply(rresults,1,function(x){diff(x[c(1,3)])/abs(x[1])}),1)))[,1:5]



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

set.seed(111)
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

cmstandard <- buildcm(classes, outputs1)
##
allscores <- apply(lut, 3, function(um){
    c(sum(cmstandard * um),
      sum(buildcm(classes, outputs1, um) * um),
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
##   standard      mixed transducer 
##  0.7571809  0.7602045  0.7629903 

pdff('../RF_transducer_gainsx', asp=1)
## tplot(x=allscores[1,1:nn2], y=allscores[2,1:nn2]-allscores[1,1:nn2], type='p', pch=16, cex=1, alpha=0.5)
tplot(x=log10(allscores[1,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, asp=1,
      xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
      xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
      yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
      ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
      xlab='utility yield from standard classification',
      ylab='utility yield from transducer & utility maximization')
abline(0,1, col=paste0(palette()[4],'88'), lwd=2, lty=1)
## tplot(x=log10(allscores[2,1:nn2]), y=log10(allscores[3,1:nn2]), type='p', pch=16, cex=0.5, alpha=0.25, col=2,
##       xticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       xlabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       yticks=log10(sort(c(1:9)*rep(10^c(-1,0),each=9))),
##       ylabels=sort(c(1:9)*rep(10^c(-1,0),each=9)),
##       xlab='utility yield from mixed method',
##       ylab='utility yield from transducer & utility maximization')
## abline(0,1, col=paste0(palette()[4],'88'), lwd=2, lty=1)
dev.off()






#########################################################
## Calculation of utility yields on altered demonstration set
## discriminative and generative mode
#########################################################
set.seed(222)
discardp <- sample(which(ddata$class==0), size=sum(ddata$class==0)-round(sum(ddata$class==1)/2))
classesb <- classes[-discardp]
outputs1b <- outputs1[-discardp]
transfoutputs1 <- cbind(X2Y[[1]](outputs1b))
colnames(transfoutputs1) <- outputcov
##
probs1b <- probs1[-discardp]
probsinv1b <- probsinv1[-discardp]
probsinv0b <- probsinv0[-discardp]

baserates <- c('0'=sum(classesb==0), '1'=sum(classesb==1))/length(classesb)
##
bayesprobs1 <- probsinv1b*baserates['1']/(probsinv0b*baserates['0'] + probsinv1b*baserates['1'])

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
      'transducer_discr'=sum(buildcm(trueclasses, probs, um) * um),
      'transducer_gener'=sum(buildcm(trueclasses, bayesprobs, um) * um)
      )
}

ulist <- list(c(1,0,0,1),
              c(1,-10,0,10),
              c(1,-100,0,100),
              c(10,0,-10,1),
              c(100,0,-100,1),
              ##
              c(1,-10,-1,10),
              c(1,-100,-1,100),
              c(10,-1,-10,1),
              c(100,-1,-100,1)
              )
##
umlist <- c(lapply(ulist, function(x){ matrix(x,2,2, byrow=T) } ),
    lapply(ulist, function(x){
        x <- x-min(x)
        x <- x/max(x)
        matrix(x,2,2, byrow=T)
    }
) )
##
ulist2 <- unlist(umlist)
dim(ulist2) <- c(4,2*length(ulist))
ulist2 <- t(ulist2)

results1 <- t(sapply(umlist, function(um){
    comparescores(trueclasses=classesb, um=um, outputs=outputs1b, probs=probs1b, bayesprobs=bayesprobs1)/length(classesb)}))
##
rresults <- round(results1,3)
##
options(width=160)
cbind(ulist2, rresults,
      'd_std'=round(apply(rresults,1,function(x){diff(x[c(1,3)])}),4),
      'd_mix'=round(apply(rresults,1,function(x){diff(x[c(2,3)])}),4),
      'rd%_std'=round(100*apply(rresults,1,function(x){diff(x[c(1,3)])/abs(x[1])}),2),
      'rd%_mix'=round(100*apply(rresults,1,function(x){diff(x[c(2,3)])/abs(x[2])}),2))

##                                         standard transducer_discr transducer_gener  d_std d_mix rd%_std rd%_mix
##  [1,]   1.000    0.000    0.000   1.000    0.835            0.918            0.959  0.124 0.041   14.85    4.47
##  [2,]   1.000    0.000  -10.000  10.000    3.746            6.609            6.849  3.103 0.240   82.84    3.63
##  [3,]   1.000    0.000 -100.000 100.000   34.483           66.080           66.667 32.184 0.587   93.33    0.89
##  [4,]  10.000  -10.000    0.000   1.000    3.797            3.333            3.840  0.043 0.507    1.13   15.21
##  [5,] 100.000 -100.000    0.000   1.000   33.428           33.333           33.458  0.030 0.125    0.09    0.38
##
##  [6,]   1.000   -1.000  -10.000  10.000    3.744            6.481            6.714  2.970 0.233   79.33    3.60
##  [7,]   1.000   -1.000 -100.000 100.000   34.481           65.266           66.333 31.852 1.067   92.38    1.63
##  [8,]  10.000  -10.000   -1.000   1.000    3.634            2.916            3.820  0.186 0.904    5.12   31.00
##  [9,] 100.000 -100.000   -1.000   1.000   33.266           32.667           33.141 -0.125 0.474   -0.38    1.45
####
## [10,]   1.000    0.000    0.000   1.000    0.835            0.918            0.959  0.124 0.041   14.85    4.47
## [11,]   0.550    0.500    0.000   1.000    0.687            0.830            0.842  0.155 0.012   22.56    1.45
## [12,]   0.505    0.500    0.000   1.000    0.672            0.830            0.833  0.161 0.003   23.96    0.36
## [13,]   1.000    0.000    0.500   0.550    0.690            0.667            0.692  0.002 0.025    0.29    3.75
## [14,]   1.000    0.000    0.500   0.505    0.667            0.667            0.667  0.000 0.000    0.00    0.00
##
## [15,]   0.550    0.450    0.000   1.000    0.687            0.824            0.836  0.149 0.012   21.69    1.46
## [16,]   0.505    0.495    0.000   1.000    0.672            0.826            0.832  0.160 0.006   23.81    0.73
## [17,]   1.000    0.000    0.450   0.550    0.682            0.646            0.691  0.009 0.045    1.32    6.97
## [18,]   1.000    0.000    0.495   0.505    0.666            0.663            0.666  0.000 0.003    0.00    0.45

t(cbind(ulist2, signif(rresults[,c(1,2,3)],3),
        round(100*apply(rresults,1,function(x){diff(x[c(1,3)])/abs(x[1])}),1),
        round(100*apply(rresults,1,function(x){diff(x[c(2,3)])/abs(x[2])}),1)
        ) )[,1:5]



cmstandard <- buildcm(classesb, outputs1b)
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
##  standard     discr     gener 
## 0.6864915 0.7401895 0.7606006 

pdff('../RF_transducer_gains_generx', asp=1)
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
pdff('../RF_transducer_gains_gener_vs_discrx', asp=1)
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

#########################################################
## 
#########################################################


#########################################################
## 
#########################################################


#########################################################
## 
#########################################################





#########################################################
## Check of Kjetil's calculations test set 2
#########################################################
kresults <- fread('RF_bayesian_prob_test2b.csv', sep=',')

X <- X2Y[[1]](data.matrix(kresults[,'pred_0']))
colnames(X) <- realCovs

ptest1 <- samplesF(Y=cbind(class=0), X=X, parmList=parmlist, inorder=F)

mptest1 <- rowMeans(ptest1)

tplot(x=kresults$prob_0, y=mptest1, type='p', xlab='Kjetil',ylab='Luca', pch='+')

discre <- (abs(mptest1-kresults[,prob_0])/mptest1)
which.max(discre)
summary(discre)*100

testhis <- thist(log10(discre));tplot(x=testhis$mids, y=testhis$density)
#########################################################

classseq <- kresults$class


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

## using Luca's probabilities
results1 <- t(sapply(umlist, function(um){
    comparescores(classseq, um=um, kresults$pred_0, kresults$pred_0, mptest1)
}))/length(classseq)
colnames(results1) <- c('standard', 'output_as_prob', 'bayes')
##
cbind(round(results1,3),t(apply(round(results1,3),1,function(x){
    mx <- max(x)
    mx <- c(x-mx)/mx*100
    names(mx) <- paste0('%diff_',names(x))
    mx
}
)))

## using Kjetil's probabilities
results2 <- t(sapply(umlist, function(um){
    comparescores(classseq, um=um, kresults$pred_0, kresults$pred_0, kresults$prob_0)
}))/length(classseq)
colnames(results2) <- c('standard', 'output_as_prob', 'bayes')
##
cbind(round(results2,3),t(apply(round(results2,3),1,function(x){
    mx <- max(x)
    mx <- c(x-mx)/mx*100
    names(mx) <- paste0('%diff_',names(x))
    mx
}
)))


kscores <- fread('scores_utility_matrix.csv', sep=',')

v1 <- results1[,1]
v2 <- results2[,1]
cbind(v1, v2, v1-v2)
max(abs(v1-v2))

v1 <- kscores$score_RF_output * length(classseq)
v2 <- results1[,1]
cbind(v1, v2, v1-v2)
max(abs(v1-v2))

v1 <- kscores$score_RF_bayesian * length(classseq)
v2 <- results1[,3]
cbind(v1, v2, v1-v2)
max(abs(v1-v2))

#########################################################
## Comparison of utility scores with & without Bayesian augmentation
## for a large number of possible utility matrices
#########################################################
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
    comparescores(classseq, um=um, kresults$pred_0, kresults$pred_0, mptest1)
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
pdff('RF_improv_umspace')
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
pdff('improvement_RF')
tplot(x=resultsl[groupp,1]/length(classseq), y=reldiff[groupp], type='p', pch=16,cex=1, col=1,
      xlab='standard method', ylab='% improvement with bayes augmentation', xlim=0:1, ylim=rg)
if(sum(groupn)>0){
    tplot(x=resultsl[groupn,1]/length(classseq), y=reldiff[groupn], type='p', pch=16,cex=1, col=2, add=T)
}
dev.off()
sum(groupp)/length(reldiff)




## tplot(x=rbind(xy[,1]), y=rbind(xy[,2]), type='p', pch=15, cex=cex, col=cols,
##       xlab='standard method', ylab='% improvement with bayes augmentation', xlim=c(-1,1), ylim=c(-1,1), asp=1)
## dev.off()

## rprec <- Inf
## reldiff <- (round(resultsl[,3],rprec)-round(resultsl[,1],rprec))/round(resultsl[,1],rprec)*100
## summary(reldiff)
## rg <- range(c(0,reldiff))
## groupp <- reldiff>= 0
## groupn <- reldiff < 0
## ##
## pdff('improvement_RF')
## tplot(x=resultsl[groupp,1]/length(classseq), y=reldiff[groupp], type='p', pch='.', col=1,
##       xlab='standard method', ylab='% improvement with bayes augmentation', xlim=0:1, ylim=rg)
## if(sum(groupn)>0){
##     tplot(x=resultsl[groupn,1]/length(classseq), y=reldiff[groupn], type='p', pch='.', col=2, add=T)
## }
## ## reldiff1 <- (round(results1[,3],rprec)-round(results1[,1],rprec))/round(results1[,1],rprec)*100
## ## tplot(x=results1[,1], y=reldiff1, type='p', pch=2, cex=2, col=3, add=T)
## dev.off()
## sum(groupp)/length(reldiff)


## groupp <- reldiff>= 0
## groupn <- reldiff < 0
## tplot(x=resultsl[,1], y=reldiff, type='p', pch=15, cex=cex, col=1,
##       xlab='standard method', ylab='% improvement with bayes augmentation', xlim=rg, ylim=rg)
## if(sum(groupn)>0){
##     tplot(x=xy[groupn,1], y=xy[groupn,2], type='p', pch=15, cex=cex, col=2, add=T)
## }


## reldiff1 <- (round(results1[,3],rprec)-round(results1[,1],rprec))/round(results1[,1],rprec)*100
## tplot(x=results1[,1], y=reldiff1, type='p', pch=2, cex=2, col=3, add=T)
## dev.off()
## sum(groupp)/length(reldiff)


## resultsl <- t(future_sapply(lum, function(um){
##     comparescores(classseq, um=um, kresults$pred_0, kresults$pred_0, mptest1)
## }))/length(classseq)
## colnames(resultsl) <- c('standard', 'output_as_prob', 'bayes')

## rprec <- Inf
## reldiff <- (round(resultsl[,3],rprec)-round(resultsl[,1],rprec))/round(resultsl[,1],rprec)*100
## summary(reldiff)
## rg <- range(c(0,reldiff))
## groupp <- reldiff>= 0
## groupn <- reldiff < 0
## ##
## pdff('improvement_RF')
## tplot(x=resultsl[groupp,1], y=reldiff[groupp], type='p', pch='.', col=1,
##       xlab='standard method', ylab='% improvement with bayes augmentation', xlim=0:1, ylim=rg)
## if(sum(groupn)>0){
##     tplot(x=resultsl[groupn,1], y=reldiff[groupn], type='p', pch='.', col=2, add=T)
## }
## ## reldiff1 <- (round(results1[,3],rprec)-round(results1[,1],rprec))/round(results1[,1],rprec)*100
## ## tplot(x=results1[,1], y=reldiff1, type='p', pch=2, cex=2, col=3, add=T)
## dev.off()
## sum(groupp)/length(reldiff)


#########################################################
## Comparison of scores for a test set with different
## proportions of classes
#########################################################
classseq <- kresults$class

props <- c(1,1)
##
csubset <- sort(c(which(classseq==1),
                  sample(which(classseq==0), round(sum(classseq==1)*props[2]/props[1]))))

probinv0 <- samplesF(X=cbind(class=0), Y=X[csubset,1, drop=F], parmList=parmlist, inorder=F)
probinv1 <- samplesF(X=cbind(class=1), Y=X[csubset,1, drop=F], parmList=parmlist, inorder=F)
##
mprobinv0 <- rowMeans(probinv0)
mprobinv1 <- rowMeans(probinv1)

probinvbase <- mprobinv0*props[1]/(mprobinv0*props[1]+mprobinv1*props[2])

## using the full MCMC results
resultsbase <- t(sapply(umlist, function(um){
    comparescores(classseq[csubset], um=um,
                  kresults$pred_0[csubset],
                  kresults$pred_0[csubset],
                  mptest1[csubset],
                  probinvbase)
}))/length(csubset)
colnames(resultsbase) <- c('standard', 'output_as_prob', 'bayes', 'inverse_bayes')

##
cbind(resultsbase[,c(1,3,4)],t(apply(resultsbase[,c(1,3,4)],1,function(x){
    mx <- max(x)
    c(x-mx)/mx*100
}
)))

resultslbase <- t(future_sapply(lum, function(um){
    comparescores(classseq[csubset], um=um, kresults$pred_0[csubset], mptest1[csubset], probinvbase)
}))/length(csubset)
colnames(resultslbase) <- c('standard', 'bayes', 'inverse_bayes')

rprec <- Inf
reldiff <- (round(resultslbase[,3],rprec)-round(resultslbase[,1],rprec))/round(resultslbase[,1],rprec)*100
summary(reldiff)
rg <- range(c(0,reldiff))
groupp <- reldiff>= 0
groupn <- reldiff < 0
##
pdff('improvement_RF_diversebalance')
tplot(x=resultsl[groupp,1], y=reldiff[groupp], type='p', pch='.', col=1,
      xlab='standard method', ylab='% improvement with inv. bayes augmentation', xlim=0:1, ylim=rg)
if(sum(groupn)>0){
    tplot(x=resultsl[groupn,1], y=reldiff[groupn], type='p', pch='.', col=2, add=T)
}
reldiff1 <- (round(resultsbase[,4],rprec)-round(resultsbase[,1],rprec))/round(resultsbase[,1],rprec)*100
tplot(x=resultsbase[,1], y=reldiff1, type='p', pch=2, cex=2, col=3, add=T)
dev.off()
sum(groupp)/length(reldiff)




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
kresults <- fread('RF_direct_prob.csv', sep=',')

X <- X2Y[[1]](data.matrix(kresults[,'preds']))
colnames(X) <- realCovs

ptest1 <- samplesF(Y=cbind(class=0), X=X, parmList=parmlist, inorder=F)

mptest1 <- rowMeans(ptest1)

discre <- (abs(mptest1-kresults[,probs_class_0])/mptest1)
which.max(discre)

summary(discre)*100
testhis <- thist(log10(discre));tplot(x=testhis$mids, y=testhis$density)
#########################################################





## Probability vs output
xgrid <- seq(0, 1, length.out=256)
ygrid <- X2Y[['prediction_lnodds']](xgrid)
##
vpoints <- cbind(prediction_lnodds=ygrid)
##
opgrid <- samplesF(Y=cbind(class=0), X=vpoints, parmList=parmlist, inorder=F)
##
qgrid <- apply(opgrid,1,function(x){quantile(x, c(1,7)/8)})
##

tplot(x=xgrid, y=cbind(rowMeans(opgrid), 1- rowMeans(opgrid)), xlab='RF % output', ylab='probability of class', ylim=c(0,1), lwd=3, family='Palatino')
legend('top', c('class 0', 'class 1'), lty=c(1,2), col=c(1,2), lwd=3, bty='n')


polygon(x=c(xgrid,rev(xgrid)), y=c(qgrid[1,],rev(qgrid[2,])), col=paste0(palette()[1],'40'), border=NA)
legend('topleft', legend=c(
                       paste0(paste0(rownames(qgrid),collapse='\u2013'), ' uncertainty')
                   ),
       lty=c('solid'), lwd=c(10),
       col=paste0(palette()[1],c('40')),
       bty='n', cex=1.25)


## Stability of probability
opgrid0 <- samplesF(X=cbind(class=0), Y=vpoints, parmList=parmlist, inorder=F)*Xjacobian[['prediction_lnodds']](xgrid)
##
qgrid0 <- apply(opgrid0,1,function(x){quantile(x, c(1,7)/8)})
##
opgrid1 <- samplesF(X=cbind(class=1), Y=vpoints, parmList=parmlist, inorder=F)*Xjacobian[['prediction_lnodds']](xgrid)
##
qgrid1 <- apply(opgrid1,1,function(x){quantile(x, c(1,7)/8)})
##
tplot(x=xgrid, y=cbind(rowMeans(opgrid0),rowMeans(opgrid1)), xlab='RF % output', ylab='probability of output', ylim=c(0,25))
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
               mu=c(shortparmlist$meanR),
               sigma=1/sqrt(c(shortparmlist$tauR))
               ), 'RF_probabilityfunction_full.csv', sep=',')


qorder <- order(c(oneparmlist$q), decreasing=FALSE)
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
## qselect <- c(oneparmlist$q)[qorder] > 2^-32
qselect <- 1:(which(cumsum(c(oneparmlist$q)[qorder]) > 1-2^-16)[1])

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
xgrid <- seq(0, 1, length.out=256)
ygrid <- X2Y[['prediction_lnodds']](xgrid)
##
vpoints <- cbind(prediction_lnodds=ygrid)
##
pgrid <- samplesF(Y=cbind(class=1), X=vpoints, parmList=shortparmlist, inorder=F)
##
## qgrid <- apply(pgrid,1,function(x){quantile(x, c(1,7)/8)})
##
tplot(x=xgrid, y=cbind(rowMeans(opgrid),pgrid), xlab='RF % output', ylab='probability of class 1', ylim=c(0,1))
## polygon(x=c(xgrid,rev(xgrid)), y=c(qgrid[1,],rev(qgrid[2,])), col=paste0(palette()[1],'40'), border=NA)
## legend('topleft', legend=c(
##                        paste0(paste0(rownames(qgrid),collapse='\u2013'), ' uncertainty')
##                    ),
##        lty=c('solid'), lwd=c(10),
##        col=paste0(palette()[1],c('40')),
##        bty='n', cex=1.25)

fwrite(data.table(w=c(shortparmlist$q),
               p=c(shortparmlist$probB),
               mu=c(shortparmlist$meanR),
               sigma=1/sqrt(c(shortparmlist$tauR))
               ), 'RF_probabilityfunction3.csv', sep=',')





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
cdata <- fread(datafile, sep=',')
alldata <- cdata[1:ndata, ..covNames]
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




orange <- c(-1,1)*ceiling(max(abs(as.matrix(cdata[,..realCovs]))))

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

