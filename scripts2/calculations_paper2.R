## Author: PGL  Porta Mana
## Created: 2022-05-01T09:38:48+0200
## Last-Updated: 2022-05-26T07:48:05+0200
################
## Calculations for the papers
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
## library('cowplot')
## library('plotly')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
library('future.apply')
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

##
## utd <- function(p,a,b,up,dup,un,dun){
##     up*a*p + un*b*(1-p) + (up-dup)*(1-a)*p + (un-dun)*(1-b)*(1-p)
## }
##
## ut <- function(p,a,b,utp,utn,ufp,ufn){
##     utp*a*p + utn*b*(1-p) + ufn*(1-a)*p + ufp*(1-b)*(1-p)
## }
##
confm1 <- function(p,a){
        matrix(c(p*a[1], p*(1-a[1]), (1-p)*(1-a[2]), (1-p)*a[2]), 2,2)
}
##
confm <- function(p,a,b){
    len <- length(p)
    if(len==1){
        matrix(c(p*a, p*(1-a), (1-p)*(1-b), (1-p)*b), 2,2)
    }else{
        aperm(array(c(p*a, p*(1-a), (1-p)*(1-b), (1-p)*b), dim=c(len,2,2)), c(2,3,1))
    }
}
##
cm2st <- function(cm){
    if(length(dim(cm))==2){
        p <- sum(cm[,1])
        c(p, cm[1,1]/p, cm[2,2]/(1-p))
    }else{
        p <- colSums(cm[,1,])
        rbind(p, cm[1,1,]/p, cm[2,2,]/(1-p))
    }
}
##
ut <- function(cm,utm){
    if(length(dim(cm))==2){
        sum(cm * utm)
    }else{
        rowSums(aperm(cm*utm, c(3,1,2)))
    }
}
##
f1score <- function(p,a,b){
    2*a*p/(2*a*p + (1-b)*(1-p) + (1-a)*p)
}
##
mcc <- function(p,a,b){
    (a*p * b*(1-p) - (1-a)*p * (1-b)*(1-p))/
        sqrt(
        (a*p + (1-b)*(1-p)) * (b*(1-p) + (1-a)*p) * p * (1-p)
        )
}
##
prec <- function(p,a,b){
    a*p/(a*p + (1-b)*(1-p))
}
##
acc <- function(p,a,b){
    a*p+b*(1-p)
}
##
bacc <- function(p,a,b){
    (a+b)/2
}
##
auc <- function(p,a,b){
    a/2 + b/2
}
##
foma <- function(p,a,b){
    sqrt(prec(p,a,b) * a)
}
##
kri <- function(p,a,b){
    raccu <- ((p*a+(1-p)*(1-b) + p)/2)^2 + (((1-p)*b+p*(1-a) + (1-p))/2)^2
    ((acc(p,a,b) + 1)/2 - raccu)/(1-raccu)
}
##
reldiff <- function(x,y){200*(x-y)/(x+y)}
##
allscores <- function(p,ab){
    out <- c(f1score(pp,ab[1],ab[2]),
      mcc(pp,ab[1],ab[2]),
      prec(pp,ab[1],ab[2]),
      acc(pp,ab[1],ab[2]),
      bacc(pp,ab[1],ab[2]),
      kri(pp,ab[1],ab[2]),
      foma(pp,ab[1],ab[2]),
      auc(pp,ab[1],ab[2]),
      ab[1],
      ab[2])
    names(out) <- c('F1', 'MCC', 'Prec', 'Acc', 'BalAcc', 'Kri', 'Fo-Ma', 'AUC', 'Rec', 'Spec')
    out
}
##
id <- diag(2)
utp <- matrix(c(1,0,0,0),2,2)
utn <- matrix(c(0,0,0,1),2,2)
ufp <- matrix(c(0,0,1,0),2,2)
ufn <- matrix(c(0,1,0,0),2,2)
##
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
##
um2xy <- function(um,norm=TRUE){
    if(norm){
        um <- um-min(um)
        um <- um/max(um)
    }
    x <- -(1-um[2,2])*(um[2,2]<1) + (1-um[1,1])*(um[1,1]<1)
    y <- -um[2,1]*(um[2,1]>0) + um[1,2]*(um[1,2]>0)
    c(x,y)
}



##########################################################################
#### Plot of wrong result of some metrics for some utility matrices
#### Version 3
##########################################################################

## xumpaper <- xumpaper-min(xumpaper)
## xumpaper <- signif(xumpaper/max(xumpaper), 2)
##
set.seed(149)
##
nn <- 3*10^3
alpha <- 0.25
#xum <- xumpaper
shape1 <- 2
shape2 <- 1
## inrange <- FALSE
## while(!inrange){
##     wxy <- um2xy(xum)+rnorm(2,0,0.1)
##     inrange <- (wxy[2] <= wxy[1]+1 & wxy[2] >= wxy[1]-1)
## }
## wum <- xy2um(wxy)
## print(wxy)
## print(wum)
##
metrlist <- list(
    'Accuracy'=acc,
    'True-positive rate'=function(p,a,b){a},
    'F1-measure'=f1score,
    'Matthews Correlation Coefficient'=mcc
    ##'Precision'=prec,
    ##'Balanced accuracy'=bacc,
    ##'Fowlkes-Mallows index'=foma,
    ##'True-negative rate'=function(p,a,b){b},
    ## ,
    ##  'Utility matrix eq. (4)'=function(p,a,b){
    ##      rowSums(aperm(confm(p,a,b)* c(xumpaper) ))
    ##  }
    )
##
lp <- rep(0.5,nn)
la1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## test <- thist(la1);tplot(x=test$breaks, y=test$density)
## summary(la1)
lb1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## la2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## lb2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
##
lcm1 <- confm(lp,la1,lb1)
## lcm2 <- confm(lp,la2,lb2)
##
## ldut <- rowSums(aperm((lcm2-lcm1)*c(xum)))
##
##
## pdff(paste0('testdiff', paste0(xum,collapse='')))
## for(i in 1:length(metrlist)){
##     metr <- metrlist[[i]]
##     ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
##     ylab <- names(metrlist)[i]
## groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
## tplot(x=list(ldut[groupok],ldut[!groupok]),
##       y=list(ldmetr[groupok],ldmetr[!groupok]), type='p', pch=c(20,17),cex=0.5,
##       alpha=0.25,
##       xlab='difference in utility', ylab=paste0('difference in ',ylab))
##     legend('topleft', bty='n', legend=paste0(
##                           'incorrectly ranked pairs: ',
##                           round(100*(1-sum(groupok)/nn)), '%'),
##            cex=1.5)
## }
## dev.off()
#################################
## apap <- 2*c(0.27,0.43)
## bpap <- 2*c(0.35, 0.32)
## lupap <- rowSums(aperm(confm(rep(0.5,2),apap,bpap)*c(xum)))
## lmpap <- metr(rep(0.5,2),apap,bpap)
## set.seed(149)
##
##
##
##
if(FALSE){
    okbs <- sapply(list(diag(2), matrix(c(1,0,0,0),2,2)), function(xum){
    print(xum)
    excl <- c(1,2)
    ldut <- rowSums(aperm((lcm1)*c(xum)))
    ##
    allmetr <- t(sapply(metrlist[-excl], function(metr){ metr(lp,la1,lb1) }))
    rgm <- apply(allmetr,1,range)
    allmetr <- (allmetr-rgm[1,])/(2*(rgm[2,]-rgm[1,]))+0.5
    rgm <- apply(allmetr,1,range)
    ldut2 <- (ldut-min(ldut))/(2*diff(range(ldut)))+0.5
    rgu <- range(ldut2)
    ##    
    system.time(test <- foreach(i=1:nn)%dopar%{
        which(apply(rbind(allmetr-allmetr[,i]<0, ldut2-ldut2[i] > 0), 2, all) |
              apply(rbind(allmetr-allmetr[,i]>0, ldut2-ldut2[i] < 0), 2, all))
    })
    ##
    tocheck <- which(sapply(test,length)>0)
    ##
    distanc <- sapply(tocheck,function(i1){
        mind <- sapply(test[[i1]],function(i2){
            min(abs(c(allmetr[,i1]-allmetr[,i2],ldut2[i1]-ldut2[i2])))
        })
        c(which.max(mind),max(mind))})
    ##
    i1 <- which.max(distanc[2,])
    i2 <- distanc[1,i1]
    c(tocheck[i1], test[[tocheck[i1]]][i2])})
}
## else{
## if(lp[1]==0.5){
##     okbs <- matrix(c(1034, 1219, 1454, 4781), 2,2)
## }else if(lp[1]==0.85){
##     okbs <- matrix(c(5899, 7508, 4781, 7411),2,2)
## }
## }
## > okbs 0.8
##      [,1] [,2]
## [1,]  138  450
## [2,] 4130 4781
##
umlist <- list(diag(2), matrix(c(1,0,0,0),2,2))
##
##
## allmetr <- t(sapply(metrlist, function(metr){ metr(lp,la1,lb1) }))
## rgm <- apply(allmetr,1,range)
## allmetr <- (allmetr-rgm[1,])/(2*(rgm[2,]-rgm[1,]))+0.5
## rgm <- apply(allmetr,1,range)
##
##     if(all(xum==diag(2))){scalingp <-c(3,3,Inf,3,1)}else{scalingp <-c(3,3,3,Inf,1)}
##     distp <- function(x){sum(abs(x)^2/scalingp )}
## okb2 <- c(which.min(apply(rbind(allmetr-rgm[2,],ldut2-rgu[1]),2, distp )),
## #         which.min(apply(rbind(allmetr-mean(rgm),ldut2-mean(rgu)),2, distp )),
##          which.min(apply(rbind(allmetr-rgm[1,],ldut2-rgu[2]),2, distp))
##          )
##
##
lem <- length(metrlist)
pdff(paste0('utility_vs_metrics2_',lp[1]),paper='a4')
pchseq <- c(5,0,1,1)
colseq <- rep(4,4)#c(3,4,6,7)
par(mfcol=c(4,2))
par(oma=c(0,0,0,0))
for(ii in 1:length(umlist)){
    xum <- umlist[[ii]]
    ldut <- rowSums(aperm((lcm1)*c(xum)))
    ldut2 <- (ldut-min(ldut))/(2*diff(range(ldut)))+0.5
    rgu <- range(ldut2)
    ##
    for(i in 1:lem){
        metr <- metrlist[[i]]
        ldmetr <- metr(lp,la1,lb1)
        ylab <- names(metrlist)[i]
        if(i==ii){#max(diff(diff(ldmetr)/diff(ldut)))<1e-10){
            ok <- c(1,2)
            pch <- NA
            col <- NA
        }else{
            ldmetr2 <- (ldmetr-min(ldmetr))/(2*diff(range(ldmetr)))+0.5
            rgm <- range(ldmetr2)
            ok <- c(which.min((ldmetr2-rgm[1])^10 + (ldut2-rgu[2])^10),
                     which.min((ldmetr2-rgm[2])^10 + (ldut2-rgu[1])^10))
            pch <- c(2,6)#rep(pchseq[i],each=2)
            col <- 2#rep(colseq[i], each=2)
        }
        ## if(i<=5){
        ##     ok1 <- ok1b
        ##     ok2 <- ok2b
        ##     pch <- 2
        ## }else{
        ##     rgm <- range(ldmetr)
        ##     ok1 <- which.min((ldmetr-rgm[1])^2 + (ldut-rgu[2])^2)
        ##     ok2 <- which.min((ldmetr-rgm[2])^2 + (ldut-rgu[1])^2)
        ##     pch <- 5
        ## }
        ## diffu <- ldut[1:(nn/2)]-ldut[(nn/2+1):nn]
        ## diffm <- ldmetr[1:(nn/2)]-ldmetr[(nn/2+1):nn]
        ## okp <- which( diffu*diffm < 0)[1]
        tplot(x=ldut, y=ldmetr,
              type='p', pch=16, cex=0.5, alpha=alpha,
              cex.axis=1.5,cex.lab=1.5, ly=3, n=5, family='Palatino',
              mar=c((if(i%/%lem!=1){4}else{6.1}),
              (if(ii%/%2!=1){4.1}else{4.1}),
              0.25,
              (if(ii==1){0.25}else{0.25})
              ),
              ylabels=(ii%/%2 != 1),
              ylab=(if(ii%/%2 != 1){ylab}else{NA}),
              lx=(i%/%lem == 1)*2.33,
              xlabels=(i%/%lem == 1),
              xlab=(if(i%/%lem == 1){
                        if(ii==1){
                            expression('yield for utility matrix '~bgroup("[",atop(1~~0,0~~1),"]"))
                        }else{
                            expression('yield for utility matrix '~bgroup("[",atop(1~~0,0~~0),"]"))
                        }
                    }else{NA})
              )
        tplot(x=rbind(ldut[ok]),
              y=rbind(ldmetr[ok]), type='p', pch=pch, cex=2, lwd=3, col=col,
              add=T)
    }
}
dev.off()


################################################################
#### Wrongly ranked pairs from metrics 
#### and from utility yields with wrong utilities
################################################################

set.seed(149)
##
nn <- 10^6
nn2 <- 3*10^3
alpha <- 0.5
shape1 <- 2
shape2 <- 1
##
metrlist <- list(
        'True-positive rate'=function(p,a,b){a},
##    'True-negative rate'=function(p,a,b){b},
    'Precision'=prec,
    'Balanced accuracy'=bacc,
    'Matthews Corr. Coeff.'=mcc,
    'Fowlkes-Mallows index'=foma,
    'F1-measure'=f1score,
    'Accuracy'=acc
)
##
lp <- runif(nn,0,1)
la1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## test <- thist(la1);tplot(x=test$breaks, y=test$density)
## summary(la1)
lb1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
la2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
lb2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
##
lcm1 <- confm(lp,la1,lb1)
lcm2 <- confm(lp,la2,lb2)
##
ipairs <- foreach(typexy=c('unif','norm'))%do%{
    print(typexy)
    if(typexy=='unif'){lxy <- runif(2*round(nn*6/3),-1,1)
    }else{ lxy <- rnorm(2*round(nn*6/3),0,1/3)}
    ##
    dim(lxy) <- c(round(nn*6/3),2)
    lxy <- lxy[lxy[,2]<lxy[,1]+1 & lxy[,2]>lxy[,1]-1 & abs(lxy[,1])<=1 & abs(lxy[,2])<=1,][1:nn,]
    ##
    lut <- apply(lxy,1,xy2um)
    dim(lut) <- c(2,2,nn)
    ##
    ldut <- rowSums(aperm((lcm2-lcm1)*lut))
    ##
    percmetrics <- sapply(metrlist, function(metr){
                ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
                groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
                round(100*(1-sum(groupok)/nn),1)
    })
    names(percmetrics) <- names(metrlist)
    ##
    errorums <- c(seq(0,50,by=10))
    percwum <- foreach(errorum=errorums, .combine=cbind)%dorng%{
        print(paste0(typexy,'_',errorum))
        wlut <- apply(lut,3,function(um){
            inrange <- FALSE
            while(!inrange){
                newum <- um+rnorm(4,0,errorum/100)
                inrange <- min(newum)>=0 & max(newum)<=1 & newum[1,1]>=newum[2,1] & newum[2,2]>=newum[1,2]
                ## newxy <- um2xy(newum,norm=T)
                ## inrange <- newxy[2]<newxy[1]+1 & newxy[2]>newxy[1]-1 & abs(newxy[1])<=1 & abs(newxy[2])<=1
            }
            newum})
        dim(wlut) <- c(2,2,nn)
        ## wlut <- apply(lxy,1,function(xy){
        ##     inrange <- FALSE
        ##     while(!inrange){
        ##         newxy <- xy+rnorm(2,0,errorum/100)
        ##         inrange <- newxy[2]<newxy[1]+1 & newxy[2]>newxy[1]-1 & abs(newxy[1])<=1 & abs(newxy[2])<=1
        ##     }
        ##     xy2um(newxy)})
        ## dim(wlut) <- c(2,2,nn)
        ##
        ldmetr <- rowSums(aperm((lcm2-lcm1)*wlut))
        ##
        groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
        print(paste0(errorum,' ',sd(wlut-lut)))
        c(sd(wlut-lut), round(100*(1-sum(groupok)/nn),1))
    }
    attr(percwum, 'rng') <- NULL
    colnames(percwum) <- errorums
    rownames(percwum) <- c('trueSD', 'percentage')
    ##
    list(metrics=percmetrics, ums=percwum)
}
names(ipairs) <- c('unif','norm')

for(typexy in c('unif','norm')){
    pdff(paste0('../increase_error2_',typexy), paper='a4r')
    datap1 <- ipairs[[typexy]]$ums
    datap2 <- ipairs[[typexy]]$metrics
    ylim <- range(c(datap1,datap2,0))
    tplot(x=datap1['trueSD',], y=datap1['percentage',], type='l', lwd=4, ylim=ylim,
          ylab='incorrectly ranked pairs/%',
          xlab='standard deviation of error in utilities',
          mar=c(4.5, 3, 0, 9)+c(1,1.1,1,1))
    abline(h=datap2, col=c(2), lwd=2, lty=2)
    mtext(text=paste0('  ',names(datap2)), side=4, at=datap2, las=1, cex=1.1, col=2)
mtext(text=paste0('(with ',(if(typexy=='unif'){'uniform'}else{'gaussian'}), ' distribution of true utility matrices)'), side=1, padj=6,cex=1.25, font=1) 
    dev.off()
}


## Plot of gaussian distribution of true UMs
tplot(x=lxy[1:nn2,1], y=lxy[1:nn2,2], type='p', pch=16, cex=0.75, col=3, alpha=0.5, xgrid=F, ygrid=F, xticks=F,yticks=F,xlab=NA,ylab=NA,asp=1, mar=rep(1,4))
tplot(x=c(-1,-1,0,1,1,0,-1), y=c(0,-1,-1,0,1,1,0), type='l', add=T, col='#000000', lwd=1)


#### Plot of examples error distribution
set.seed(149)
errorum <- 11
lxy <- rep(c(0.5,0.5),each=nn2)
dim(lxy) <- c(nn2,2)
lut <- apply(lxy,1,xy2um)
dim(lut) <- c(2,2,nn2)
wlut <- apply(lut,3,function(um){
    inrange <- FALSE
    while(!inrange){
        newum <- um+rnorm(4,0,errorum/100)
        inrange <- min(newum)>=0 & max(newum)<=1 & newum[1,1]>=newum[2,1] & newum[2,2]>=newum[1,2]
        ## newxy <- um2xy(newum,norm=T)
        ## inrange <- newxy[2]<newxy[1]+1 & newxy[2]>newxy[1]-1 & abs(newxy[1])<=1 & abs(newxy[2])<=1
    }
    newum})
dim(wlut) <- c(2,2,nn2)
sd(wlut-lut)
wlxy <- t(apply(wlut,3,um2xy))
##
tplot(x=wlxy[1:1000,1], y=wlxy[1:1000,2], type='p', pch=17, cex=0.75, col=1, alpha=0.5, xgrid=F, ygrid=F, xticks=F,yticks=F,xlab=NA,ylab=NA,asp=1, mar=rep(1,4), xlim=c(-1,1), ylim=c(-1,1))
tplot(x=lxy[1,1], y=lxy[1,2], type='p', add=T, col='#000000', cex=3, pch=18)
##
errorum <- 50
lxy <- rep(c(-0.5,-0.5),each=nn2)
dim(lxy) <- c(nn2,2)
lut <- apply(lxy,1,xy2um)
dim(lut) <- c(2,2,nn2)
wlut <- apply(lut,3,function(um){
    inrange <- FALSE
    while(!inrange){
        newum <- um+rnorm(4,0,errorum/100)
        inrange <- min(newum)>=0 & max(newum)<=1 & newum[1,1]>=newum[2,1] & newum[2,2]>=newum[1,2]
        ## newxy <- um2xy(newum,norm=T)
        ## inrange <- newxy[2]<newxy[1]+1 & newxy[2]>newxy[1]-1 & abs(newxy[1])<=1 & abs(newxy[2])<=1
    }
    newum})
dim(wlut) <- c(2,2,nn2)
sd(wlut-lut)
wlxy <- t(apply(wlut,3,um2xy))
##
tplot(x=wlxy[1:1000,1], y=wlxy[1:1000,2], type='p', pch=15, cex=0.75, col=2, alpha=0.5, xgrid=F, ygrid=F, xticks=F,yticks=F,xlab=NA,ylab=NA,asp=1, mar=rep(1,4), xlim=c(-1,1), ylim=c(-1,1), add=T)
tplot(x=lxy[1,1], y=lxy[1,2], type='p', add=T, col='#000000', cex=3, pch=18)
##
tplot(x=c(-1,-1,0,1,1,0,-1), y=c(0,-1,-1,0,1,1,0), type='l', add=T, col='#000000', lwd=1)



######################################################################
#### 
######################################################################



######################################################################
#### 
######################################################################



######################################################################
#### 
######################################################################



######################################################################
#### 
######################################################################



######################################################################
#### 
######################################################################



######################################################################
#### 
######################################################################



######################################################################
#### 
######################################################################




######################################################################
#### Example of bad sampling from dataset, need of inverse probability
######################################################################

mm <- 80
aa <- 1:(mm-1)
testx <- outer(aa, aa, function(y,x){x})
testx2 <- outer(aa, aa, function(y,x){mm-x})
testy <- outer(aa, aa, function(y,x){y})
testy2 <- outer(aa, aa, function(y,x){mm-y})
sele <- which(
    testx*testy/mm==round(testx*testy/mm) &
    testx2*testy/mm==round(testx2*testy/mm) &
    testx*testy2/mm==round(testx*testy2/mm) &
    testx2*testy2/mm==round(testx2*testy2/mm) &
    mm*testx*testy/(testx*testy+testx2*testy2)==round(mm*testx*testy/(testx*testy+testx2*testy2)) &
    mm*testx*testy2/(testx*testy2+testx2*testy)==round(mm*testx*testy2/(testx*testy2+testx2*testy)) &
   ,
    arr.ind=T
)
sele <- sele[sele[,1]!=mm/2 & sele[,2]!=mm/2,]
sele

cho <- 5
y <- sele[cho,1]
x <- sele[cho,2]
##
m1 <- matrix(c(aa[x], mm-aa[x], rep(mm/4,2), mm-aa[x], aa[x]),2,3)
#m1 <- matrix(c(aa[x], mm-aa[x], mm-aa[x], aa[x]),2,2)
m2 <- m1*c(aa[y], mm-aa[y])
m1
m2
m1/rowSums(m1)
m2/rowSums(m2)
t(t(m1)/colSums(m1))
t(t(m2)/colSums(m2))
##      [,1] [,2] [,3]
## [1,]   24   20   56
## [2,]   56   20   24
## >      [,1] [,2] [,3]
## [1,]  480  400 1120
## [2,] 3360 1200 1440
## >      [,1] [,2] [,3]
## [1,] 0.24  0.2 0.56
## [2,] 0.56  0.2 0.24
## >      [,1] [,2] [,3]
## [1,] 0.24  0.2 0.56
## [2,] 0.56  0.2 0.24
## >      [,1] [,2] [,3]
## [1,]  0.3  0.5  0.7
## [2,]  0.7  0.5  0.3
## >       [,1] [,2]   [,3]
## [1,] 0.125 0.25 0.4375
## [2,] 0.875 0.75 0.5625



####################################################################


lo <- round((1-0.54)/0.02) + 1
qq <- seq(0.54,1,length.out=lo)
listab <- cbind('Rec'=rep(qq,lo), 'Spec'=rep(qq,each=lo))
pp <- 0.5
##
listcm <- future_apply(listab, 1, function(ab){
    confm1(pp,ab)
})
dim(listcm) <- c(2,2,lo*lo)

##
listscores <- t(future_apply(listab, 1, function(ab){
    c(f1score(pp,ab[1],ab[2]),
      mcc(pp,ab[1],ab[2]),
      prec(pp,ab[1],ab[2]),
      acc(pp,ab[1],ab[2]),
      bacc(pp,ab[1],ab[2]))
}))
colnames(listscores) <- c('F1', 'MCC', 'Prec', 'Acc', 'BalAcc')

lo2 <- round((0.9+0.9)/0.1) + 1
ss <- round(seq(-0.9,0.9,length.out=lo2),2)
xy <- cbind(rep(ss,lo2), rep(ss,each=lo2))
xy <- xy[xy[,2]<xy[,1]+1 & xy[,2]>xy[,1]-1,]
lxy <- nrow(xy)
##tplot(x=xy[,1], y=xy[,2], type='p')
##
listum <- future_apply(xy, 1, xy2um)
dim(listum) <- c(2,2,lxy)


listutil <- future_apply(cbind(rep(1:lxy,lo*lo), rep(1:(lo*lo),each=lxy)), 1,
                       function(ab){
                           sum(listum[,,ab[1]] * listcm[,,ab[2]])
                       })
dim(listutil) <- c(lxy,lo*lo)


diffscores <- future_sapply(
    1:ncol(listscores), function(i){
        200*c( outer(listscores[,i], listscores[,i], '-')/
               outer(listscores[,i], listscores[,i], '+') )
    })
colnames(diffscores) <- colnames(listscores)
##
diffscores2 <- future_sapply(
    1:ncol(listab), function(i){
        200*c( outer(listab[,i], listab[,i], '-')/
               outer(listab[,i], listab[,i], '+') )
    })
colnames(diffscores2) <- colnames(listab)


diffutilities <- future_sapply(
    1:lxy, function(i){
        200*c(outer(listutil[i,], listutil[i,], '-')/
              outer(listutil[i,], listutil[i,], '+') )
    })

dsselectp <- apply(diffscores,1,function(x){all(x>0)})
dsselectp2 <- sapply(1:nrow(diffscores),function(x){all(abs(diffscores[x,])>abs(diffscores2[x,'Spec']))})

selectn <- future_apply(diffutilities, 2, function(x){
    (x<0 & dsselectp & diffscores2[,'Rec']>0 & dsselectp2) 
})

dim(selectn) <- c(lo*lo,lo*lo, lxy)
selectn <- which(selectn, arr.ind=T)
colnames(selectn) <- c('icm1','icm2','ium')

## subsel <- which((listum[1,1,selectn[,'ium']]==listum[2,1,selectn[,'ium']] |
##                           listum[1,2,selectn[,'ium']]==listum[2,2,selectn[,'ium']])
##                   ##       &
##                   ##       !(listcm[1,1,selectn[,2]]==listcm[2,1,selectn[,2]] |
##                   ## listcm[1,2,selectn[,2]]==listcm[2,2,selectn[,2]]) &
##                   ##       !(listcm[1,1,selectn[,1]]==listcm[2,1,selectn[,1]] |
##                   ##         listcm[1,2,selectn[,1]]==listcm[2,2,selectn[,1]])
##                   )

subsel <- selectn[which( future_apply(selectn,1,function(it){
    (min(diag(listum[,,it['ium']])) > 1*sum(listcm[,,it['icm1']]*listum[,,it['ium']])) &
        abs(reldiff(sum(listcm[,,it['icm1']]*listum[,,it['ium']]), sum(listcm[,,it['icm2']]*listum[,,it['ium']])))>=2
})),]
nrow(subsel)

ok <- subsel[which.min(abs(diffscores2[(subsel[,'icm1']-1)*lo^2+subsel[,'icm2'], 'Spec'])),]
ok1 <- (ok['icm1']-1)*lo^2 + ok['icm2']
##
tcm1 <- listcm[,,ok['icm1']]
print('CM1')
tcm1
tcm2 <- listcm[,,ok['icm2']]
print('CM2')
tcm2
tab1 <- listab[ok['icm1'],]
tab2 <- listab[ok['icm2'],]
##
tum <- listum[,,ok['ium']]
tut1 <- sum(tum*tcm1)
tut2 <- sum(tum*tcm2)
##
tum <- (tum-mean(c(tut1,tut2)))
tut1 <- sum(tum*tcm1)
tut2 <- sum(tum*tcm2)
tum <- tum*(10^ceiling(-log10(tut2)))
##
tum <- round((tum-min(tum))/2) - 335 # factory
## tum <- round((tum-min(tum))/2) - 335 + c(25,-25,0,0)# factory inverse
##tum <- round((tum-min(tum))/2) # medical
print('UM')
tum
tut1 <- sum(tum*tcm1)
print('utility CM1')
signif(tut1/c(1,12),4)
tut2 <- sum(tum*tcm2)
print('utility CM2')
signif(tut2/c(1,12),4)
print('% Dutility CM1 - CM2')
signif(reldiff(tut1,tut2),2)
##
tsc1 <- allscores(pp, listab[ok['icm1'],])
print('scores CM1')
signif(tsc1,2)
tsc2 <- allscores(pp, listab[ok['icm2'],])
print('scores CM2')
signif(tsc2,2)
print('% Dscores CM1 - CM2')
signif(reldiff(tsc1,tsc2),2)

signif(listscores[ok['icm1'],],2)
signif(listscores[ok['icm2'],],2)
signif(diffscores[ok1,],2)

#### Paper example 1
## [1] "CM1"
## >      [,1] [,2]
## [1,] 0.43 0.18
## [2,] 0.07 0.32
## > > [1] "CM2"
## >      [,1] [,2]
## [1,] 0.27 0.15
## [2,] 0.23 0.35
## > > > > > > > > > > > > > > > > [1] "UM"
## >      [,1] [,2]
## [1,]   15 -335
## [2,]  -35  165
## > > [1] "utility CM1"
## > [1] -3.5000 -0.2917
## > > [1] "utility CM2"
## > [1] 3.5000 0.2917
## > [1] "% Dutility CM1 - CM2"
## > [1] 2.6e+17
## > > > [1] "scores CM1"
## >   F1    MCC   Prec    Acc BalAcc    Kri  Fo-Ma    AUC    Rec   Spec 
##   0.77   0.51   0.70   0.75   0.75   0.75   0.78   0.75   0.86   0.64 
## > > [1] "scores CM2"
## >   F1    MCC   Prec    Acc BalAcc    Kri  Fo-Ma    AUC    Rec   Spec 
##   0.59   0.24   0.64   0.62   0.62   0.62   0.59   0.62   0.54   0.70 
## > [1] "% Dscores CM1 - CM2"
## >     F1    MCC   Prec    Acc BalAcc    Kri  Fo-Ma    AUC    Rec   Spec 
##   28.0   71.0    9.2   19.0   19.0   19.0   28.0   19.0   46.0   -9.0 
##
#### Paper example 1 small change
## [1] "CM1"
## >      [,1] [,2]
## [1,] 0.43 0.18
## [2,] 0.07 0.32
## > > [1] "CM2"
## >      [,1] [,2]
## [1,] 0.27 0.15
## [2,] 0.23 0.35
## > > > > > > > > > > > > > > > > [1] "UM"
## >      [,1] [,2]
## [1,]   45 -335
## [2,]  -65  165
## > > [1] "utility CM1"
## > [1] 7.3000 0.6083
## > > [1] "utility CM2"
## > [1] 4.7000 0.3917
## > [1] "% Dutility CM1 - CM2"
## > [1] 43
## > > > [1] "scores CM1"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   0.77   0.51   0.70   0.75   0.75   0.75   0.75   0.86   0.64 
## > > [1] "scores CM2"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   0.59   0.24   0.64   0.62   0.62   0.62   0.62   0.54   0.70 
## > [1] "% Dscores CM1 - CM2"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   28.0   71.0    9.2   19.0   19.0   19.0   19.0   46.0   -9.0 
##
#### Paper example 1 small change II
## [1] "CM1"
## >      [,1] [,2]
## [1,] 0.43 0.18
## [2,] 0.07 0.32
## > > [1] "CM2"
## >      [,1] [,2]
## [1,] 0.27 0.15
## [2,] 0.23 0.35
## > > > > > > > > > > > > > > > > [1] "UM"
## >      [,1] [,2]
## [1,]   40 -335
## [2,]  -60  165
## > > [1] "utility CM1"
## > [1] 5.5000 0.4583
## > > [1] "utility CM2"
## > [1] 4.500 0.375
## > [1] "% Dutility CM1 - CM2"
## > [1] 20
## > > > [1] "scores CM1"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   0.77   0.51   0.70   0.75   0.75   0.75   0.75   0.86   0.64 
## > > [1] "scores CM2"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   0.59   0.24   0.64   0.62   0.62   0.62   0.62   0.54   0.70 
## > [1] "% Dscores CM1 - CM2"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   28.0   71.0    9.2   19.0   19.0   19.0   19.0   46.0   -9.0 

#### Paper example 1 medicine
## [1] "CM1"
## >      [,1] [,2]
## [1,] 0.43 0.18
## [2,] 0.07 0.32
## > > [1] "CM2"
## >      [,1] [,2]
## [1,] 0.27 0.15
## [2,] 0.23 0.35
## > > > > > > > > > > > > > > > > [1] "UM"
## >      [,1] [,2]
## [1,]  350    0
## [2,]  300  500
## > > [1] "utility CM1"
## > [1] 331.50  27.62
## > > [1] "utility CM2"
## > [1] 338.50  28.21
## > [1] "% Dutility CM1 - CM2"
## > [1] -2.1
## > > > [1] "scores CM1"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   0.77   0.51   0.70   0.75   0.75   0.75   0.75   0.86   0.64 
## > > [1] "scores CM2"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   0.59   0.24   0.64   0.62   0.62   0.62   0.62   0.54   0.70 
## > [1] "% Dscores CM1 - CM2"
## >     F1    MCC   Prec    Acc BalAcc    Kri    AUC    Rec   Spec 
##   28.0   71.0    9.2   19.0   19.0   19.0   19.0   46.0   -9.0 




##########################################################################
#### Plot of wrong result of some metrics for some utility matrices
##########################################################################

xumpaper <- matrix(c(15,-35,-335,165), 2,2)
## xumpaper <- xumpaper-min(xumpaper)
## xumpaper <- signif(xumpaper/max(xumpaper), 2)
##
set.seed(149)
##
nn <- 10^4
xum <- diag(2)
#xum <- xumpaper
shape1 <- 2
shape2 <- 1
## inrange <- FALSE
## while(!inrange){
##     wxy <- um2xy(xum)+rnorm(2,0,0.1)
##     inrange <- (wxy[2] <= wxy[1]+1 & wxy[2] >= wxy[1]-1)
## }
## wum <- xy2um(wxy)
## print(wxy)
## print(wum)
##
metrlist <- list('F1-measure'=f1score,
              'MCC'=mcc,
              'Precision'=prec,
              'Balanced accuracy'=bacc,
              'Fowlkes-Mallows index'=foma,
              'True-positive rate'=function(p,a,b){a},
              'True-negative rate'=function(p,a,b){b},
              'Accuracy'=acc,
              'Utility matrix eq. (4)'=function(p,a,b){
                  rowSums(aperm(confm(p,a,b)* c(xumpaper) ))
              }
              )
##
lp <- runif(nn,0,1)
la1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## test <- thist(la1);tplot(x=test$breaks, y=test$density)
## summary(la1)
lb1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## la2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## lb2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
##
lcm1 <- confm(lp,la1,lb1)
## lcm2 <- confm(lp,la2,lb2)
##
## ldut <- rowSums(aperm((lcm2-lcm1)*c(xum)))
##
##
## pdff(paste0('testdiff', paste0(xum,collapse='')))
## for(i in 1:length(metrlist)){
##     metr <- metrlist[[i]]
##     ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
##     ylab <- names(metrlist)[i]
## groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
## tplot(x=list(ldut[groupok],ldut[!groupok]),
##       y=list(ldmetr[groupok],ldmetr[!groupok]), type='p', pch=c(20,17),cex=0.5,
##       alpha=0.25,
##       xlab='difference in utility', ylab=paste0('difference in ',ylab))
##     legend('topleft', bty='n', legend=paste0(
##                           'incorrectly ranked pairs: ',
##                           round(100*(1-sum(groupok)/nn)), '%'),
##            cex=1.5)
## }
## dev.off()
#################################
## apap <- 2*c(0.27,0.43)
## bpap <- 2*c(0.35, 0.32)
## lupap <- rowSums(aperm(confm(rep(0.5,2),apap,bpap)*c(xum)))
## lmpap <- metr(rep(0.5,2),apap,bpap)
## set.seed(149)
##
##
##

okbs <- sapply(list(diag(2), xumpaper, matrix(c(1,0,0,0),2,2)), function(xum){
    print(xum)
    excl <- (if(all(xum==diag(2))){8}else{9})
    ldut <- rowSums(aperm((lcm1)*c(xum)))
    ##
    allmetr <- t(sapply(metrlist[-excl], function(metr){ metr(lp,la1,lb1) }))
    rgm <- apply(allmetr,1,range)
    allmetr <- (allmetr-rgm[1,])/(2*(rgm[2,]-rgm[1,]))+0.5
    rgm <- apply(allmetr,1,range)
    ldut2 <- (ldut-min(ldut))/(2*diff(range(ldut)))+0.5
    rgu <- range(ldut2)
    ##    
    test <- sapply(1:nn, function(i){
        which(apply(rbind(allmetr-allmetr[,i]<0, ldut2-ldut2[i] > 0), 2, all) |
              apply(rbind(allmetr-allmetr[,i]>0, ldut2-ldut2[i] < 0), 2, all))
    })
    ##
    tocheck <- which(sapply(test,length)>0)
    ##
    distanc <- sapply(tocheck,function(i1){
        mind <- sapply(test[[i1]],function(i2){
            min(abs(c(allmetr[,i1]-allmetr[,i2],ldut2[i1]-ldut2[i2])))
        })
        c(which.max(mind),max(mind))})
    ##
    i1 <- which.max(distanc[2,])
    i2 <- distanc[1,i1]
    c(tocheck[i1], test[[tocheck[i1]]][i2])})

## okbs <- matrix(c(3754, 7902, 4751, 9865), 2,2)


umlist <- list(diag(2), xumpaper, matrix(c(1,0,0,0),2,2))
##
for(ii in 1:length(umlist)){
    xum <- umlist[[ii]]
ldut <- rowSums(aperm((lcm1)*c(xum)))
##
allmetr <- t(sapply(metrlist, function(metr){ metr(lp,la1,lb1) }))
rgm <- apply(allmetr,1,range)
allmetr <- (allmetr-rgm[1,])/(2*(rgm[2,]-rgm[1,]))+0.5
rgm <- apply(allmetr,1,range)
ldut2 <- (ldut-min(ldut))/(2*diff(range(ldut)))+0.5
rgu <- range(ldut2)
okb <- okbs[,ii]
    ## okb <- c(which.min(colSums(((allmetr-rgm[2,])^2)/c(9,9,9,9,9,9,9,9,5)) + (ldut2-rgu[1])^2),
##          which.min(colSums(((allmetr-rgm[1,])^2)/c(9,9,9,9,9,9,9,9,5)) + (ldut2-rgu[2])^2)
##          )
## okb <- c(which.min(apply(rbind(allmetr-rgm[2,],ldut2-rgu[1]),2,
##                          function(x){sum(abs(x)^2/c(9,9,9,9,9,9,9,9,9,1))})),
##          which.min(apply(rbind(allmetr-rgm[1,],ldut2-rgu[2]),2,
##                          function(x){sum(abs(x)^2/c(9,9,9,9,9,9,9,9,9,1))})))
## tryn <- 6
## okb <- c(which(!is.na(test))[tryn],test[which(!is.na(test))[tryn]])
##
##
pdff(paste0('utility_vs_metrics_', paste0(xum,collapse='')),paper='a4r')
##pngf(paste0('utility_vs_metrics_', paste0(xum,collapse='')))
exc <- 0
pchseq <- c(5,0,1,6)
colseq <- c(4,6,7)
par(mfrow=c(3,3))
par(oma=c(0,0,0,0))
for(i in 1:length(metrlist)){
    metr <- metrlist[[i]]
    ldmetr <- metr(lp,la1,lb1)
    ylab <- names(metrlist)[i]
    if( diff(ldut[okb])*diff(ldmetr[okb]) < 0){
        ok <- okb
        pch <- 2
        col <- 2
    }else if(max(diff(diff(ldmetr)/diff(ldut)))==0){
        ok <- okb
        pch <- NA
        col <- NA
        }else{
        ldmetr2 <- (ldmetr-min(ldmetr))/(2*diff(range(ldmetr)))+0.5
        rgm <- range(ldmetr2)
        ok <- c(which.min((ldmetr2-rgm[1])^2 + (ldut2-rgu[2])^2),
                which.min((ldmetr2-rgm[2])^2 + (ldut2-rgu[1])^2))
        exc <- exc + 1
        pch <- pchseq[exc]
        col <- colseq[exc]
    }
    ## if(i<=5){
    ##     ok1 <- ok1b
    ##     ok2 <- ok2b
    ##     pch <- 2
    ## }else{
    ##     rgm <- range(ldmetr)
    ##     ok1 <- which.min((ldmetr-rgm[1])^2 + (ldut-rgu[2])^2)
    ##     ok2 <- which.min((ldmetr-rgm[2])^2 + (ldut-rgu[1])^2)
    ##     pch <- 5
    ## }
    ## diffu <- ldut[1:(nn/2)]-ldut[(nn/2+1):nn]
    ## diffm <- ldmetr[1:(nn/2)]-ldmetr[(nn/2+1):nn]
    ## okp <- which( diffu*diffm < 0)[1]
tplot(x=ldut,
      y=ldmetr, type='p', pch=20, cex=0.5, alpha=0.66, cex.axis=1.5,cex.lab=1.5, ly=3.5, n=5, family='Palatino',
      mar=(if(i<6){c(1, 5.2, 2, 1)}else{c(4.1, 5.2, 2, 1)}),
      xlabels=(i>6),
      xlab=(if(i>6){'utility'}else{NA}),
      ylab=ylab)
tplot(x=rbind(ldut[ok]),
      y=rbind(ldmetr[ok]), type='p', pch=pch,cex=2, lwd=3, col=col,
      add=T)
#legend('topleft', legend=sum(groupok)/nn, bty='n')
}
dev.off()
}

##########################################################################
#### Plot of wrong result of some metrics for some utility matrices
#### Version 2
##########################################################################

## xumpaper <- xumpaper-min(xumpaper)
## xumpaper <- signif(xumpaper/max(xumpaper), 2)
##
set.seed(149)
##
nn <- 10^4
#xum <- xumpaper
shape1 <- 2
shape2 <- 1
## inrange <- FALSE
## while(!inrange){
##     wxy <- um2xy(xum)+rnorm(2,0,0.1)
##     inrange <- (wxy[2] <= wxy[1]+1 & wxy[2] >= wxy[1]-1)
## }
## wum <- xy2um(wxy)
## print(wxy)
## print(wum)
##
metrlist <- list(
    'Accuracy'=acc,
    'True-positive rate'=function(p,a,b){a},
    'F1-measure'=f1score,
    'Matthews Correlation Coefficient'=mcc
    ##'Precision'=prec,
    ##'Balanced accuracy'=bacc,
    ##'Fowlkes-Mallows index'=foma,
    ##'True-negative rate'=function(p,a,b){b},
    ## ,
    ##  'Utility matrix eq. (4)'=function(p,a,b){
    ##      rowSums(aperm(confm(p,a,b)* c(xumpaper) ))
    ##  }
    )
##
lp <- rep(0.8,nn)
la1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## test <- thist(la1);tplot(x=test$breaks, y=test$density)
## summary(la1)
lb1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## la2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## lb2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
##
lcm1 <- confm(lp,la1,lb1)
## lcm2 <- confm(lp,la2,lb2)
##
## ldut <- rowSums(aperm((lcm2-lcm1)*c(xum)))
##
##
## pdff(paste0('testdiff', paste0(xum,collapse='')))
## for(i in 1:length(metrlist)){
##     metr <- metrlist[[i]]
##     ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
##     ylab <- names(metrlist)[i]
## groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
## tplot(x=list(ldut[groupok],ldut[!groupok]),
##       y=list(ldmetr[groupok],ldmetr[!groupok]), type='p', pch=c(20,17),cex=0.5,
##       alpha=0.25,
##       xlab='difference in utility', ylab=paste0('difference in ',ylab))
##     legend('topleft', bty='n', legend=paste0(
##                           'incorrectly ranked pairs: ',
##                           round(100*(1-sum(groupok)/nn)), '%'),
##            cex=1.5)
## }
## dev.off()
#################################
## apap <- 2*c(0.27,0.43)
## bpap <- 2*c(0.35, 0.32)
## lupap <- rowSums(aperm(confm(rep(0.5,2),apap,bpap)*c(xum)))
## lmpap <- metr(rep(0.5,2),apap,bpap)
## set.seed(149)
##
##
##

if(FALSE){
    okbs <- sapply(list(diag(2), matrix(c(1,0,0,0),2,2)), function(xum){
    print(xum)
    excl <- (if(all(xum==diag(2))){1}else{2})
    ldut <- rowSums(aperm((lcm1)*c(xum)))
    ##
    allmetr <- t(sapply(metrlist[-excl], function(metr){ metr(lp,la1,lb1) }))
    rgm <- apply(allmetr,1,range)
    allmetr <- (allmetr-rgm[1,])/(2*(rgm[2,]-rgm[1,]))+0.5
    rgm <- apply(allmetr,1,range)
    ldut2 <- (ldut-min(ldut))/(2*diff(range(ldut)))+0.5
    rgu <- range(ldut2)
    ##    
    system.time(test <- foreach(i=1:nn)%dopar%{
        which(apply(rbind(allmetr-allmetr[,i]<0, ldut2-ldut2[i] > 0), 2, all) |
              apply(rbind(allmetr-allmetr[,i]>0, ldut2-ldut2[i] < 0), 2, all))
    })
    ##
    tocheck <- which(sapply(test,length)>0)
    ##
    distanc <- sapply(tocheck,function(i1){
        mind <- sapply(test[[i1]],function(i2){
            min(abs(c(allmetr[,i1]-allmetr[,i2],ldut2[i1]-ldut2[i2])))
        })
        c(which.max(mind),max(mind))})
    ##
    i1 <- which.max(distanc[2,])
    i2 <- distanc[1,i1]
    c(tocheck[i1], test[[tocheck[i1]]][i2])})
}else{
if(lp[1]==0.5){
    okbs <- matrix(c(1034, 1219, 1454, 4781), 2,2)
}else if(lp[1]==0.8){
    okbs <- matrix(c(138, 4130, 450, 4781),2,2)
}
}

umlist <- list(diag(2), matrix(c(1,0,0,0),2,2))
##
##
## allmetr <- t(sapply(metrlist, function(metr){ metr(lp,la1,lb1) }))
## rgm <- apply(allmetr,1,range)
## allmetr <- (allmetr-rgm[1,])/(2*(rgm[2,]-rgm[1,]))+0.5
## rgm <- apply(allmetr,1,range)
##
##     if(all(xum==diag(2))){scalingp <-c(3,3,Inf,3,1)}else{scalingp <-c(3,3,3,Inf,1)}
##     distp <- function(x){sum(abs(x)^2/scalingp )}
## okb2 <- c(which.min(apply(rbind(allmetr-rgm[2,],ldut2-rgu[1]),2, distp )),
## #         which.min(apply(rbind(allmetr-mean(rgm),ldut2-mean(rgu)),2, distp )),
##          which.min(apply(rbind(allmetr-rgm[1,],ldut2-rgu[2]),2, distp))
##          )
##
##
lem <- length(metrlist)
pdff(paste0('utility_vs_metrics_',lp[1]),paper='a4')
pchseq <- c(5,0,1,1)
colseq <- rep(4,4)#c(3,4,6,7)
par(mfcol=c(4,2))
par(oma=c(0,0,0,0))
for(ii in 1:length(umlist)){
    xum <- umlist[[ii]]
    ldut <- rowSums(aperm((lcm1)*c(xum)))
    ldut2 <- (ldut-min(ldut))/(2*diff(range(ldut)))+0.5
    rgu <- range(ldut2)
    ##
    for(i in 1:lem){
        metr <- metrlist[[i]]
        ldmetr <- metr(lp,la1,lb1)
        ylab <- names(metrlist)[i]
        if(max(diff(diff(ldmetr)/diff(ldut)))<1e-11){
            ok1 <- ok2 <- c(1,2)
            pch <- NA
            col <- NA
        }else{
            ok1 <- okbs[,ii]
            ldmetr2 <- (ldmetr-min(ldmetr))/(2*diff(range(ldmetr)))+0.5
            rgm <- range(ldmetr2)
            ok2 <- c(which.min((ldmetr2-rgm[1])^10 + (ldut2-rgu[2])^10),
                     which.min((ldmetr2-rgm[2])^10 + (ldut2-rgu[1])^10))
            pch <- rep(c(2+4*(ii%/%2==1),pchseq[i]),each=2)
            col <- rep(c(2,colseq[i]), each=2)
        }
        ## if(i<=5){
        ##     ok1 <- ok1b
        ##     ok2 <- ok2b
        ##     pch <- 2
        ## }else{
        ##     rgm <- range(ldmetr)
        ##     ok1 <- which.min((ldmetr-rgm[1])^2 + (ldut-rgu[2])^2)
        ##     ok2 <- which.min((ldmetr-rgm[2])^2 + (ldut-rgu[1])^2)
        ##     pch <- 5
        ## }
        ## diffu <- ldut[1:(nn/2)]-ldut[(nn/2+1):nn]
        ## diffm <- ldmetr[1:(nn/2)]-ldmetr[(nn/2+1):nn]
        ## okp <- which( diffu*diffm < 0)[1]
        tplot(x=ldut, y=ldmetr,
              type='p', pch=20, cex=0.5, alpha=0.66,
              cex.axis=1.5,cex.lab=1.5, ly=3, n=5, family='Palatino',
              mar=c((if(i%/%lem!=1){4}else{6.1}),
              (if(ii%/%2!=1){4.1}else{4.1}),
              0.25,
              (if(ii==1){0.25}else{0.25})
              ),
              ylabels=(ii%/%2 != 1),
              ylab=(if(ii%/%2 != 1){ylab}else{NA}),
              lx=(i%/%lem == 1)*2.33,
              xlabels=(i%/%lem == 1),
              xlab=(if(i%/%lem == 1){
                        if(ii==1){
                            expression('yield for utility matrix '~bgroup("[",atop(1~~0,0~~1),"]"))
                        }else{
                            expression('yield for utility matrix '~bgroup("[",atop(1~~0,0~~0),"]"))
                        }
                    }else{NA})
              )
        tplot(x=rbind(ldut[c(ok1,ok2)]),
              y=rbind(ldmetr[c(ok1,ok2)]), type='p', pch=pch, cex=2, lwd=3, col=col,
              add=T)
    }
}
dev.off()


################################################################
#### Plot of metric values for pairs of confusion matrices and
#### draws of utility matrices
#### Version 2
################################################################

set.seed(149)
##
nn <- 10^6
nn2 <- 3*10^3
alpha <- 0.5
shape1 <- 2
shape2 <- 1
##
metrlist <- list(
##        'True-positive rate'=function(p,a,b){a},
##    'True-negative rate'=function(p,a,b){b},
    'Precision'=prec,
    'Balanced accuracy'=bacc,
    'Matthews Corr. Coeff.'=mcc,
    'Fowlkes-Mallows index'=foma,
    'F1-measure'=f1score,
    'Accuracy'=acc
)
##
lp <- runif(nn,0,1)
la1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## test <- thist(la1);tplot(x=test$breaks, y=test$density)
## summary(la1)
lb1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
la2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
lb2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
##
lcm1 <- confm(lp,la1,lb1)
lcm2 <- confm(lp,la2,lb2)
##
lxy <- runif(2*round(nn*6/3),-1,1)
dim(lxy) <- c(round(nn*6/3),2)
lxy <- lxy[lxy[,2]<lxy[,1]+1 & lxy[,2]>lxy[,1]-1,][1:nn,]
##
lut <- apply(lxy,1,xy2um)
dim(lut) <- c(2,2,nn)
##
ldut <- rowSums(aperm((lcm2-lcm1)*lut))
##
func <- identity#function(x){plogis(x*10)*2-1}
##
lem <- length(metrlist)
errorums <- c(20,10)
pdff(paste0('incorrectscores2_',errorums[1]), paper='a4')
par(mfrow=c(lem/2+1,2))
par(oma=c(0,0,0,0))
for(i in 1:(lem+2)){
    if(i <= lem){
        metr <- metrlist[[i]]
        ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
        ylab <- names(metrlist)[i]
    } else {
        errorum <- errorums[i-lem]
        wlxy <- t(apply(lxy,1, function(xy){
            xy2 <- xy+rnorm(2,0,errorum/100)
            tempum <- xy2um(xy2)
            tempum <- tempum-min(tempum)
            tempum <- tempum/max(tempum)
            um2xy(tempum)
        }))
##
        wlut <- apply(wlxy,1,xy2um)
        dim(wlut) <- c(2,2,nn)
        ##
        wldut <- rowSums(aperm((lcm2-lcm1)*wlut))
        ##
        summary(sapply(1:dim(lut)[3],function(i){mean(abs(lut[,,i]-wlut[,,i]))/mean(lut[,,i])}))
        ldmetr <- wldut
        ylab <- paste0('utility, ',errorum,'% incorrect utilities')
    }
    groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
    groupok2 <- (ldut[1:nn2]>0 & ldmetr[1:nn2]>0) | (ldut[1:nn2]<0 & ldmetr[1:nn2]<0)
    tplot(x=list(func(ldut[1:nn2][groupok2]),func(ldut[1:nn2][!groupok2])),
          y=list(func(ldmetr[1:nn2][groupok2]),func(ldmetr[1:nn2][!groupok2])),
          family='Palatino', 
          type='p', pch=c(20,17), cex=c(0.5,0.65),
          mar=(if(i<=lem){c(3, 5.2, 2, 2)}else{c(4.1, 5.2, 2, 2)}),
          xlabels=(i>lem), 
          xlab=(if(i>lem){bquote(Delta~'utility yield')}else{NA}),
          ##col.lab=c('black',(if(i>lem){3}else{'black'})),
          col=c(1,(if(i>lem){4}else{2})),
          alpha=alpha,
          cex.axis=1.5, ly=3.5, n=5, cex.lab=1.5,
          ylab=bquote(Delta~.(ylab)))
    text(x=min(ldut[1:nn2]),y=max(ldmetr[1:nn2]), xpd=T, offset=0,
         labels=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn),1), '%'),
         adj=c(0.05,-0.5), cex=1.5, col=(if(i>lem){4}else{2})
         )
    ## legend(x=min(ldut),y=max(ldmetr), bty='n', text.col=2, cex=1.5, xjust=0,
    ##        legend=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn)), '%'))
}
dev.off()

################################################################
#### Plot of metric values for pairs of confusion matrices and
#### draws of utility matrices
#### Version 3
################################################################

set.seed(149)
##
nn <- 10^6
nn2 <- 3*10^3
alpha <- 0.5
shape1 <- 2
shape2 <- 1
##
metrlist <- list(
##        'True-positive rate'=function(p,a,b){a},
##    'True-negative rate'=function(p,a,b){b},
    'Precision'=prec,
    'Balanced accuracy'=bacc,
    'Matthews Corr. Coeff.'=mcc,
    'Fowlkes-Mallows index'=foma,
    'F1-measure'=f1score,
    'Accuracy'=acc
)
##
lp <- runif(nn,0,1)
la1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## test <- thist(la1);tplot(x=test$breaks, y=test$density)
## summary(la1)
lb1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
la2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
lb2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
##
lcm1 <- confm(lp,la1,lb1)
lcm2 <- confm(lp,la2,lb2)
##
lxy <- rnorm(2*round(nn*6/3),0,1/2)
dim(lxy) <- c(round(nn*6/3),2)
lxy <- lxy[lxy[,2]<lxy[,1]+1 & lxy[,2]>lxy[,1]-1,][1:nn,]
##
lut <- apply(lxy,1,xy2um)
dim(lut) <- c(2,2,nn)
##
ldut <- rowSums(aperm((lcm2-lcm1)*lut))
##
func <- identity#function(x){plogis(x*10)*2-1}
##
lem <- length(metrlist)
errorums <- c(25,10)
pdff(paste0('incorrectscores2norm2_',errorums[1]), paper='a4')
par(mfrow=c(lem/2+1,2))
par(oma=c(0,0,0,0))
for(i in 1:(lem+2)){
    if(i <= lem){
        metr <- metrlist[[i]]
        ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
        ylab <- names(metrlist)[i]
    } else {
        errorum <- errorums[i-lem]
        wlxy <- t(apply(lxy,1, function(xy){
            xy2 <- xy+rnorm(2,0,errorum/100)
            tempum <- xy2um(xy2)
            tempum <- tempum-min(tempum)
            tempum <- tempum/max(tempum)
            um2xy(tempum)
        }))
##
        wlut <- apply(wlxy,1,xy2um)
        dim(wlut) <- c(2,2,nn)
        ##
        wldut <- rowSums(aperm((lcm2-lcm1)*wlut))
        ##
        summary(sapply(1:dim(lut)[3],function(i){mean(abs(lut[,,i]-wlut[,,i]))/mean(lut[,,i])}))
        ldmetr <- wldut
        ylab <- paste0('utility, ',errorum,'% incorrect utilities')
    }
    groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
    groupok2 <- (ldut[1:nn2]>0 & ldmetr[1:nn2]>0) | (ldut[1:nn2]<0 & ldmetr[1:nn2]<0)
    tplot(x=list(func(ldut[1:nn2][groupok2]),func(ldut[1:nn2][!groupok2])),
          y=list(func(ldmetr[1:nn2][groupok2]),func(ldmetr[1:nn2][!groupok2])),
          family='Palatino', 
          type='p', pch=c(20,17), cex=c(0.5,0.65),
          mar=(if(i<=lem){c(3, 5.2, 2, 2)}else{c(4.1, 5.2, 2, 2)}),
          xlabels=(i>lem), 
          xlab=(if(i>lem){bquote(Delta~'utility yield')}else{NA}),
          ##col.lab=c('black',(if(i>lem){3}else{'black'})),
          col=c(1,(if(i>lem){4}else{2})),
          alpha=alpha,
          cex.axis=1.5, ly=3.5, n=5, cex.lab=1.5,
          ylab=bquote(Delta~.(ylab)))
    text(x=min(ldut[1:nn2]),y=max(ldmetr[1:nn2]), xpd=T, offset=0,
         labels=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn),1), '%'),
         adj=c(0.05,-0.5), cex=1.5, col=(if(i>lem){4}else{2})
         )
    ## legend(x=min(ldut),y=max(ldmetr), bty='n', text.col=2, cex=1.5, xjust=0,
    ##        legend=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn)), '%'))
}
dev.off()


################################################################
#### Plot of metric values for pairs of confusion matrices and
#### draws of utility matrices
#### Version 4
################################################################

set.seed(149)
##
nn <- 10^6
nn2 <- 3*10^3
alpha <- 0.5
shape1 <- 2
shape2 <- 1
##
metrlist <- list(
        'True-positive rate'=function(p,a,b){a},
##    'True-negative rate'=function(p,a,b){b},
    'Precision'=prec,
    'Balanced accuracy'=bacc,
    'Matthews Corr. Coeff.'=mcc,
    'Fowlkes-Mallows index'=foma,
    'F1-measure'=f1score,
    'Accuracy'=acc
)
##
lp <- runif(nn,0,1)
la1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## test <- thist(la1);tplot(x=test$breaks, y=test$density)
## summary(la1)
lb1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
la2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
lb2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
##
lcm1 <- confm(lp,la1,lb1)
lcm2 <- confm(lp,la2,lb2)
##
for(typexy in c('unif','norm')){
    if(typexy=='unif'){lxy <- runif(2*round(nn*6/3),-1,1)
    }else{ lxy <- rnorm(2*round(nn*6/3),0,1/3)}
dim(lxy) <- c(round(nn*6/3),2)
lxy <- lxy[lxy[,2]<lxy[,1]+1 & lxy[,2]>lxy[,1]-1,][1:nn,]
##
lut <- apply(lxy,1,xy2um)
dim(lut) <- c(2,2,nn)
##
ldut <- rowSums(aperm((lcm2-lcm1)*lut))
##
func <- identity#function(x){plogis(x*10)*2-1}
##
lem <- length(metrlist)
errorums <- c(11)
pdff(paste0('../incorrectscores4',typexy,'-',paste0(errorums,collapse='_')), paper='a4')
par(mfrow=c(lem/2+1,2))
par(oma=c(0,0,0,0))
for(i in 1:(lem+length(errorums))){
    if(i <= lem){
        metr <- metrlist[[i]]
        ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
        ylab <- names(metrlist)[i]
    } else {
        errorum <- errorums[i-lem]
        wlut <- apply(lut,3,function(um){
            inrange <- FALSE
            while(!inrange){
                newum <- um+rnorm(4,0,errorum/100)
                inrange <- min(newum)>=0 & max(newum)<=1 & newum[1,1]>=newum[2,1] & newum[2,2]>=newum[1,2]
                ## newxy <- um2xy(newum,norm=T)
                ## inrange <- newxy[2]<newxy[1]+1 & newxy[2]>newxy[1]-1 & abs(newxy[1])<=1 & abs(newxy[2])<=1
            }
            newum})
        dim(wlut) <- c(2,2,nn)
        print(paste0(typexy,'_',sd(wlut-lut)))
        ## wlut <- apply(lxy,1,function(xy){
        ##     xy <- lxy[i,]
        ##     inrange <- FALSE
        ##     while(!inrange){
        ##         newxy <- xy+rnorm(2,0,errorum/100)
        ##         inrange <- newxy[2]<newxy[1]+1 & newxy[2]>newxy[1]-1 & abs(newxy[1])<=1 & abs(newxy[2])<=1
        ##     }
        ##     xy2um(newxy)})
        ## dim(wlut) <- c(2,2,nn)
        ## wlut <- lut+rnorm(length(lut),0,errorum/100)
        ##
        wldut <- rowSums(aperm((lcm2-lcm1)*wlut))
        ##
        ## summary(sapply(1:dim(lut)[3],function(i){mean(abs(lut[,,i]-wlut[,,i]))/mean(lut[,,i])}))
        ldmetr <- wldut
        ylab <- paste0('utility, ',round(sd(wlut-lut),1),'-std error')
    }
    groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
    groupok2 <- (ldut[1:nn2]>0 & ldmetr[1:nn2]>0) | (ldut[1:nn2]<0 & ldmetr[1:nn2]<0)
    tplot(x=list(func(ldut[1:nn2][groupok2]),func(ldut[1:nn2][!groupok2])),
          y=list(func(ldmetr[1:nn2][groupok2]),func(ldmetr[1:nn2][!groupok2])),
          family='Palatino', 
          type='p', pch=c(16,17), cex=c(0.5,0.65),
          mar=(if(i<lem+length(errorums)-2){c(3, 5.2, 2, 2)}else{c(4.1, 5.2, 2, 2)}),
          xlabels=(i>lem+length(errorums)-2), 
          xlab=(if(i>(lem+length(errorums)-2)){bquote(Delta~'utility yield')}else{NA}),
          ##col.lab=c('black',(if(i>lem){3}else{'black'})),
          col=c(1,(if(i>lem){4}else{2})),
          alpha=alpha,
          cex.axis=1.5, ly=3.5, n=5, cex.lab=1.5,
          ylab=bquote(Delta~.(ylab)))
    text(x=min(ldut[1:nn2]),y=max(ldmetr[1:nn2]), xpd=T, offset=0,
         labels=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn),1), '%'),
         adj=c(0.05,-0.5), cex=1.5, col=(if(i>lem){4}else{2})
         )
    ## legend(x=min(ldut),y=max(ldmetr), bty='n', text.col=2, cex=1.5, xjust=0,
    ##        legend=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn)), '%'))
}
dev.off()
}




################################################################
#### AUC example
################################################################

mm <- 3
tplot(x=cbind(c(0, 0.1, 1), c(0,0.5,1)),
      y=cbind(c(0, 0.1*mm, 1), c(0,1,1)),
      asp=1, xgrid=F, ygrid=F, ylab=expression(t), xlab=expression(f),
      xlim=c(0,1), ylim=c(0,1),
      lwd=3)
for(qq in seq(-2,2,0.5)){
    abline(a=qq, b=rep(mm,3), col=7, ylim=c(0,1))
}
tplot(x=c(0,1,1,0,0), y=c(0,0,1,1,0), type='l', add=T, col='#000000', lwd=0.5)


mm <- 1/4
inte1 <- 0.78
inte2 <- 1-mm
dinte <- inte1-inte2
##
line1 <- xspline(x=c(0, 0.7, 1), y=c(0, 0.7*mm+0.8, 1), shape=c(0,0.5,0), lwd=3, border=palette()[1], draw=F)
line2 <- xspline(x=c(0,0.1,1), y=c(0,0.75,1), shape=c(0,0.5,0), lwd=3, lty=2, border=palette()[2], draw=F)
##
tplot(x=list(line1$x, line2$x), y=list(line1$y, line2$y),
    ## x=cbind(c(0, 0.8, 1), c(0,0.1,1)),
    ##   y=cbind(c(0, 0.8*mm+0.75, 1), c(0,0.6,1)),
      asp=1, xgrid=F, ygrid=F,
      ylab=bquote('true-positive rate'~~italic(t)),
      xlab=bquote('false-positive rate'~~italic(f)),
      xlim=c(0,1), ylim=c(0,1),
      mar=c(3.5, 4, 0, 0)+c(1,1.1,1,0),
      lwd=3, family='Palatino')
cols <- c(7,7,2,1,7,7)
ltys <- c(7,7,6,5,7,7)
lwds <- c(1,1,3,3,1,1)/2
k <- 0
for(qq in seq(inte2-2*dinte,inte1+2*dinte,by=dinte)){
    k <- k+1
    abline(a=qq, b=rep(mm,3), col=cols[k], lwd=lwds[k])
}
tplot(x=c(0,1,1,0,0), y=c(0,0,1,1,0), type='l', add=T, col='#000000', lwd=0.5)

## (y-1)/(0.8*mm+0.75 - 1) = (x-1)/(0.8-1)
## y= 

##     0.6/(mm*0.1) =k
## 1 = mm +k
    
##     (y-1)/(0.6 - 1) = (x-1)/(0.1-1)
## y= - (0.6 - 1)/(0.1-1) + 1






################################################################
#### Plot of metric values for pairs of confusion matrices and
#### draws of utility matrices
################################################################


set.seed(149)
##
nn <- 10^5
nn2 <- 3*10^3
alpha <- 0.25
shape1 <- 2
shape2 <- 1
##
metrlist <- list('F1-measure'=f1score,
              'Matthews Corr. Coeff.'=mcc,
              'Precision'=prec,
              'Balanced accuracy'=bacc,
              'Fowlkes-Mallows index'=foma,
              'Accuracy'=acc,
              'True-positive rate'=function(p,a,b){a},
              'True-negative rate'=function(p,a,b){b})
##
lp <- runif(nn,0,1)
la1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
## test <- thist(la1);tplot(x=test$breaks, y=test$density)
## summary(la1)
lb1 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
la2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
lb2 <- 0.5+0.5*rbeta(nn, shape1=shape1, shape2=shape2)
##
lcm1 <- confm(lp,la1,lb1)
lcm2 <- confm(lp,la2,lb2)
##
for(errorum in c(10,15,20,25)){
    print(errorum)
lxy <- runif(2*round(nn*6/3),-1,1)
dim(lxy) <- c(round(nn*6/3),2)
lxy <- lxy[lxy[,2]<lxy[,1]+1 & lxy[,2]>lxy[,1]-1,][1:nn,]
##
## wlxy <- t(apply(lxy,1, function(xy){
##     inrange <- FALSE
##     k <- 0
##     while((!inrange) & k<100){
##         k <- k+1
##         xy2 <- xy+rnorm(2,0,errorum/100)
##         inrange <- (xy2[2] <= xy2[1]+1 & xy2[2] >= xy2[1]-1)
##     }
##     c(xy2,k)}))
## if(max(wlxy[,3])>=99){print('WARNING')}
## wlxy <- wlxy[,1:2]
wlxy <- t(apply(lxy,1, function(xy){
    xy2 <- xy+rnorm(2,0,errorum/100)
    tempum <- xy2um(xy2)
    tempum <- tempum-min(tempum)
    tempum <- tempum/max(tempum)
    um2xy(tempum)
    }))
##
lut <- apply(lxy,1,xy2um)
dim(lut) <- c(2,2,nn)
##
ldut <- rowSums(aperm((lcm2-lcm1)*lut))
##
wlut <- apply(wlxy,1,xy2um)
dim(wlut) <- c(2,2,nn)
##
wldut <- rowSums(aperm((lcm2-lcm1)*wlut))
##
summary(sapply(1:dim(lut)[3],function(i){mean(abs(lut[,,i]-wlut[,,i]))/mean(lut[,,i])}))
##
func <- identity#function(x){plogis(x*10)*2-1}
##
    endi0 <- 1
    endi <- length(metrlist)+1
    pdff(paste0('incorrectscores_',errorum), paper='a4r')
    par(mfrow=c(3,3))
    par(oma=c(0,0,0,0))
    for(i in endi0:endi){
        if(i < length(metrlist)+1){
            metr <- metrlist[[i]]
    ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
    ylab <- names(metrlist)[i]
    } else {
    ldmetr <- wldut
    ylab <- paste0('utility, ',errorum,'% incorrect utilities')
    }
        groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
        groupok2 <- (ldut[1:nn2]>0 & ldmetr[1:nn2]>0) | (ldut[1:nn2]<0 & ldmetr[1:nn2]<0)
        tplot(x=list(func(ldut[1:nn2][groupok2]),func(ldut[1:nn2][!groupok2])),
              family='Palatino',
              y=list(func(ldmetr[1:nn2][groupok2]),func(ldmetr[1:nn2][!groupok2])),
              type='p', pch=c(20,17),cex=c(0.5,0.65),
              mar=(if(i<6){c(1, 5.2, 2, 1)}else{c(4.1, 5.2, 2, 1)}),
              xlabels=(i>6),
##              xlab=(if(i>6){'diff. in utility'}else{NA}),
              xlab=(if(i>6){bquote(Delta~'utility yield')}else{NA}),
              alpha=alpha, cex.axis=1.5, cex.lab=1.5, ly=3.5, n=5,
              ylab=bquote(Delta~.(ylab)))
        text(x=min(ldut[1:nn2]),y=max(ldmetr[1:nn2]), xpd=T, offset=0,
             labels=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn)), '%'),
             adj=c(0.05,-0.5), cex=1.5, col=2
             )
        ## legend(x=min(ldut),y=max(ldmetr), bty='n', text.col=2, cex=1.5, xjust=0,
        ##        legend=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn)), '%'))
    }
    dev.off()
}











for(j in 1:3){
    endi0 <- (j-1)*3+1
    endi <- j*3
    pdff(paste0('draws_',j), paper='a4', width=11.7, height=16.5)
    par(mfrow=c(3,1))
    par(oma=c(0,0,0,0))
    for(i in endi0:endi){
        if(i < length(metrlist)+1){
            metr <- metrlist[[i]]
    ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
    ylab <- names(metrlist)[i]
    } else {
    ldmetr <- wldut
    ylab <- paste0('utility, ',errorum,'% incorrect utilities')
    }
        groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
        tplot(x=list(func(ldut[1:nn2][groupok2]),func(ldut[1:nn2][!groupok2])), 
              y=list(func(ldmetr[1:nn2][groupok2]),func(ldmetr[1:nn2][!groupok2])), type='p', pch=c(20,17),cex=c(0.5,0.65),
              alpha=0.66, xlabels=(!(i<endi)),
              mar=(if(i<endi){c(1, 5.2, 2, 0)}else{c(4.1, 5.2, 2, 0)}),
              xlab=(if(i<endi){NA}else{'diff. in utility'}),
              ylab=paste0('diff. in ',ylab))
        legend('topleft', bty='n', text.col=2, cex=1.5,
               legend=paste0('incorrectly ranked pairs: ', round(100*(1-sum(groupok)/nn)), '%'))
    }
    dev.off()
}


## pdff(paste0('draws'))
## for(i in 1:(length(metrlist)+1)){
##     if(i < length(metrlist)+1){
##     metr <- metrlist[[i]]
##     ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
##     ylab <- names(metrlist)[i]
##     } else {
##     ldmetr <- wldut
##     ylab <- paste0('utility using ',errorum,'% incorrect utilities')
##     }
##     groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
##     tplot(x=list(ldut[groupok],ldut[!groupok]),
##       y=list(ldmetr[groupok],ldmetr[!groupok]), type='p', pch=c(20,17),cex=0.5,
##       alpha=0.5,
##       xlab='difference in utility', ylab=paste0('difference in ',ylab))
##     legend('topleft', bty='n', legend=paste0(
##                           'incorrectly ranked pairs: ',
##                           round(100*(1-sum(groupok)/nn)), '%'),
##            cex=1.5)
## }
## dev.off()

## pdff(paste0('draws2'), paper='a4', width=11.7, height=16.5)
## endi <- 5
## par(mfrow=c(endi,1))
## par(oma=c(0,0,0,0))
## for(i in 1:endi){
##     if(i < length(metrlist)+1){
##     metr <- metrlist[[i]]
##     ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
##     ylab <- names(metrlist)[i]
##     } else {
##     ldmetr <- wldut
##     ylab <- paste0('utility using ',errorum,'% incorrect utilities')
##     }
##     groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
##     tplot(x=list(ldut[groupok],ldut[!groupok]), 
##       y=list(ldmetr[groupok],ldmetr[!groupok]), type='p', pch=c(20,17),cex=0.5,
##       alpha=0.5, xlabels=(!(i<endi)),
##       mar=(if(i<endi){c(1, 5.2, 1, 0)}else{c(4.1, 5.2, 1, 0)}),
##       xlab=(if(i<endi){NA}else{'difference in utility'}),
##       ylab=paste0('difference in ',ylab))
##     legend('topleft', bty='n', legend=paste0(
##                           'incorrectly ranked pairs: ',
##                           round(100*(1-sum(groupok)/nn)), '%'),
##            cex=1.5)
## }
## dev.off()
## ##
## pdff(paste0('draws2b'), paper='a4', width=11.7, height=16.5)
## endi0 <- endi+1
## endi <- length(metrlist)+1
## par(mfrow=c(endi-endi0+1,1))
## par(oma=c(0,0,0,0))
## for(i in endi0:endi){
##     if(i < length(metrlist)+1){
##     metr <- metrlist[[i]]
##     ldmetr <- metr(lp,la2,lb2)-metr(lp,la1,lb1)
##     ylab <- names(metrlist)[i]
##     } else {
##     ldmetr <- wldut
##     ylab <- paste0('utility using ',errorum,'% incorrect utilities')
##     }
##     groupok <- (ldut>0 & ldmetr>0) | (ldut<0 & ldmetr<0)
##     tplot(x=list(ldut[groupok],ldut[!groupok]), 
##       y=list(ldmetr[groupok],ldmetr[!groupok]), type='p', pch=c(20,17),cex=0.5,
##       alpha=0.5, xlabels=(!(i<endi)),
##       mar=(if(i<endi){c(1, 5.2, 1, 0)}else{c(4.1, 5.2, 1, 0)}),
##       xlab=(if(i<endi){NA}else{'difference in utility'}),
##       ylab=paste0('difference in ',ylab))
##     legend('topleft', bty='n', legend=paste0(
##                           'incorrectly ranked pairs: ',
##                           round(100*(1-sum(groupok)/nn)), '%'),
##            cex=1.5)
## }
## dev.off()





#################################
apap <- 2*c(0.27,0.43)
bpap <- 2*c(0.35, 0.32)
lupap <- rowSums(aperm(confm(rep(0.5,2),apap,bpap)*c(xum)))
lmpap <- metr(rep(0.5,2),apap,bpap)
set.seed(149)
##
##
##
ldut <- rowSums(aperm((lcm1)*c(xum)))
rgu <- range(ldut)
allmetr <- sapply(metrlist, function(metr){ metr(lp,la1,lb1) })
rgm <- apply(allmetr,2,range)
okb <- c(which.min(colSums((t(allmetr)-rgm[2,])^2)/ncol(allmetr) + (ldut-rgu[1])^2),
         which.min(colSums((t(allmetr)-rgm[1,])^2)/ncol(allmetr) + (ldut-rgu[2])^2)
         )




##
##
pdff(paste0('testsingle', paste0(xum,collapse='')))
exc <- 0
pchseq <- c(5,0,1,6)
for(i in 1:length(metrlist)){
    metr <- metrlist[[i]]
    ldmetr <- metr(lp,la1,lb1)
    ylab <- names(metrlist)[i]
    if( diff(ldut[okb])*diff(ldmetr[okb]) < 0){
        ok <- okb
        pch <- 2
        col <- 2
    }else{
        rgm <- range(ldmetr)
        ok <- c(which.min((ldmetr-rgm[1])^2 + (ldut-rgu[2])^2),
                which.min((ldmetr-rgm[2])^2 + (ldut-rgu[1])^2))
        exc <- exc + 1
        pch <- pchseq[exc]
        col <- 2 + exc
    }
    ## if(i<=5){
    ##     ok1 <- ok1b
    ##     ok2 <- ok2b
    ##     pch <- 2
    ## }else{
    ##     rgm <- range(ldmetr)
    ##     ok1 <- which.min((ldmetr-rgm[1])^2 + (ldut-rgu[2])^2)
    ##     ok2 <- which.min((ldmetr-rgm[2])^2 + (ldut-rgu[1])^2)
    ##     pch <- 5
    ## }
    ## diffu <- ldut[1:(nn/2)]-ldut[(nn/2+1):nn]
    ## diffm <- ldmetr[1:(nn/2)]-ldmetr[(nn/2+1):nn]
    ## okp <- which( diffu*diffm < 0)[1]
tplot(x=ldut,
      y=ldmetr, type='p', pch=20, cex=0.5, alpha=0.25,
      xlab='utility', ylab=ylab)
tplot(x=rbind(ldut[ok]),
      y=rbind(ldmetr[ok]), type='p', pch=pch,cex=3, lwd=5, col=col,
      add=T)
#legend('topleft', legend=sum(groupok)/nn, bty='n')
}
dev.off()


####################################################################
####################################################################
####################################################################





select <- future_apply(diffutilities, 2, function(x){
    (x>0 & dsselectn) | (x<0 & dsselectp)
})

dim(select) <- c(lo*lo,lo*lo, lxy)



rgs <- range(diffscores)
rgu <- range(sapply(1:lxy, function(umn){
    c(outer(utilities[umn,], utilities[umn,], '-'))
}))
for(umn in 1: lxy){
diffutilities <- c(outer(utilities[umn,], utilities[umn,], '-'))
##
agree <- (diffscores>0 & diffutilities>0) | (diffscores<0 & diffutilities<0) | (diffscores==0 & diffutilities==0)
tplot(x=diffutilities[agree], y=diffscores[agree], type='p', pch=16,cex=0.5, alpha=7/8, col=1, xlim=rgu, ylim=rgs, add=(umn>1))
tplot(x=diffutilities[!agree], y=diffscores[!agree], type='p', pch='.', col=2, add=T)
}



####################################################################

nn <- 64
qqa <- runif(nn*2, 0.5,1)
qqb <- runif(nn*2, 0.5,1)
pp <- 0.5
##
listscores <- t(future_apply(cbind(qqa,qqb), 1, function(ab){
    c(f1score(pp,ab[1],ab[2]),
      mcc(pp,ab[1],ab[2]),
      prec(pp,ab[1],ab[2]),
      acc(pp,ab[1],ab[2]),
      bacc(pp,ab[1],ab[2]))
}))
colnames(listscores) <- c('F1', 'MCC', 'Prec', 'Acc', 'BalAcc')

##
listcm <- future_apply(cbind(qqa,qqb), 1, function(ab){
    confm(pp,ab[1],ab[2])
})
dim(listcm) <- c(2,2,nn*2)

lou <- round(128*8/6)
ssa <- runif(lou,-1,1)
ssb <- runif(lou,-1,1)
xy <- cbind(ssa,ssb)
xy <- xy[xy[,2]<=xy[,1]+1 & xy[,2]>=xy[,1]-1,]
lxy <- nrow(xy)
##tplot(x=xy[,1], y=xy[,2], type='p')
##
id <- diag(2)
utp <- matrix(c(1,0,0,0),2,2)
utn <- matrix(c(0,0,0,1),2,2)
ufp <- matrix(c(0,0,1,0),2,2)
ufn <- matrix(c(0,1,0,0),2,2)
##
listum <- future_apply(xy[1:lxy,], 1, function(ab){
        id + (ab[1]<0)*ab[1]*utn - (ab[1]>0)*ab[1]*utp +
        (ab[2]>0)*ab[2]*ufp - (ab[2]<0)*ab[2]*ufn
})
dim(listum) <- c(2,2,lxy)


utilities <- future_apply(cbind(rep(1:lxy,nn*2), rep(1:(nn*2),each=lxy)), 1,
                       function(ab){
                           sum(listum[,,ab[1]] * listcm[,,ab[2]])
                       })
dim(utilities) <- c(lxy,nn*2)

metr <- 2
diffscores <- c(outer(listscores[1:nn,metr], listscores[(nn+1):(nn*2),metr], '-'))
rgs <- range(diffscores)
rgu <- range(sapply(1:lxy, function(umn){
    c(outer(utilities[umn,1:nn], utilities[umn,(nn+1):(nn*2)], '-'))
}))
ini <- FALSE
sumagree <- 0
for(umn in 1:lxy){
diffutilities <- c(outer(utilities[umn,1:nn], utilities[umn,(nn+1):(nn*2)], '-'))
##
agree <- (diffscores>0 & diffutilities>0) | (diffscores<0 & diffutilities<0) | (diffscores==0 & diffutilities==0)
tplot(x=diffutilities[agree], y=diffscores[agree], type='p', pch='.',, col=1, xlim=rgu, ylim=rgs, add=ini)
tplot(x=diffutilities[!agree], y=diffscores[!agree], type='p', pch='.', col=2, add=T)
ini <- TRUE
sumagree <- sumagree+sum(agree)
}
sumagree/(lxy*nn)




future_apply(xy, 1, function(coo){
    outer(qq, qq, FUN=function(alpha,beta){
        
    })
})



test <- future_sapply(1:(lo*lo*lxy), function(i){
    coo <- xy[(ceiling(i - 1)%%lxy)+1,]
    alpha <- qq[(ceiling(i/lxy - 1)%%lo)+1]
    beta <- qq[(ceiling(i/(lxy*lo) - 1)%%lo)+1]
    ##
    um <- id + (coo[1]<0)*coo[1]*utn - (coo[1]>0)*coo[1]*utp +
        (coo[2]>0)*coo[2]*ufp - (coo[2]<0)*coo[2]*ufn
    ##
    ut <- sum(um * matrix(c(pp*alpha, pp*(1-alpha), (1-pp)*(1-beta), (1-pp)*beta),2,2))


}


resultsl <- t(future_apply(xy, 1, function(coo){
    if(coo[2]<=coo[1]+1 & coo[2]>=coo[1]-1){
    um <- id + (coo[1]<0)*coo[1]*utn - (coo[1]>0)*coo[1]*utp +
        (coo[2]>0)*coo[2]*ufp - (coo[2]<0)*coo[2]*ufn
    ##
    comparescores(classseq, um=um, kresults$sigmoid, kresults$sigmoid, mptest1)
    }else{rep(NA,3)}
}))
colnames(resultsl) <- c('standard', 'output_as_prob', 'bayes')


set.seed(149)
nn <- 10^4
lp <- runif(nn, 0.5, 1)
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
lup <- runif(nn, 0, 1)
ldup <- runif(nn, 0, 1)
lun <- runif(nn, 0, 1)
ldun <- runif(nn, 0, 1)
##

ut1 <- ut(lp,la1,lb1,lup,ldup,lun,ldun)
ut2 <- ut(lp,la2,lb2,lup,ldup,lun,ldun)
Du <- ut2-ut1
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
Dac <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)

tplot(x=Du, y=Df,
      type='p',pch='.', xlab='U(1)-U(2)', ylab='F1(1)-F1(2)')

tplot(x=Du, y=Dm,
      type='p',pch='.', xlab='U(1)-U(2)', ylab='MCC(1)-MCC(2)')


## sele <- which( (ut2 > 0 & ut1 < 0 & Df < 0 & Dm < 0) | (ut2 < 0 & ut1 > 0 & Df > 0 & Dm > 0) )

## sele <- which( (ut2 > 0 & ut1 < 0 & Df < 0 & Dm < 0 & Dpr < 0 & Dac < 0) | (ut2 < 0 & ut1 > 0 & Df > 0 & Dm > 0 & Dpr > 0 & Dac > 0) )
sele <- which(Du < 0 & Df > 0 & Dm > 0 & Dpr > 0 & Dac > 0)
length(sele)

tplot(x=Du[sele], y=Df[sele],
      type='p',pch='.', xlab='U(1)-U(2)', ylab='F1(1)-F1(2)')

tplot(x=Du[sele], y=Dm[sele],
      type='p',pch='.', xlab='U(1)-U(2)', ylab='MCC(1)-MCC(2)')

## okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele]) +ldup[sele]+ldun[sele])]
## okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele])+abs(Dpr[sele])+abs(Dac[sele]) +ldup[sele]+ldun[sele]+lp[sele])]
okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele])+abs(Dpr[sele])+abs(Dac[sele]) +ldup[sele]+ldun[sele]-lp[sele])]
##
print('P+:')
lp[okpoint]
c(la1[okpoint], lb1[okpoint])
c(la2[okpoint], lb2[okpoint])
##
print('CM1:')
matrix(c(la1[okpoint]*lp[okpoint],(1-la1[okpoint])*lp[okpoint],
         (1-lb1[okpoint])*(1-lp[okpoint]),lb1[okpoint]*(1-lp[okpoint])), 2,2)
print('CM2:')
matrix(c(la2[okpoint]*lp[okpoint],(1-la2[okpoint])*lp[okpoint],
         (1-lb2[okpoint])*(1-lp[okpoint]),lb2[okpoint]*(1-lp[okpoint])), 2,2)
print('UM:')
matrix(c(lup[okpoint], lup[okpoint]-ldup[okpoint],
         lun[okpoint]-ldun[okpoint], lun[okpoint]),2,2)
##
print('utilities:')
c(ut(lp,la1,lb1,lup,ldup,lun,ldun)[okpoint], ut(lp,la2,lb2,lup,ldup,lun,ldun)[okpoint])
print('F1 scores:')
c(f1score(lp,la1,lb1)[okpoint], f1score(lp,la2,lb2)[okpoint])
print('MCCs:')
c(mcc(lp,la1,lb1)[okpoint], mcc(lp,la2,lb2)[okpoint])
print('Precs:')
c(prec(lp,la1,lb1)[okpoint], prec(lp,la2,lb2)[okpoint])
print('Accs:')
c(acc(lp,la1,lb1)[okpoint], acc(lp,la2,lb2)[okpoint])

tp <- 0.55
ta1 <- 0.6
tb1 <- 0.6
ta2 <- 0.99
tb2 <- 0.51
tup <- 5
tdup <- tup - 4
tun <- 0
tdun <- tun -10
##
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
matrix(c(ta1*tp,(1-ta1)*tp,
         (1-tb1)*(1-tp),tb1*(1-tp)), 2,2)
print('CM2:')
matrix(c(ta2*tp,(1-ta2)*tp,
         (1-tb2)*(1-tp),tb2*(1-tp)), 2,2)
print('UM:')
matrix(c(tup, tup-tdup,
         tun-tdun, tun),2,2)
##
print('utilities:')
c(ut(tp,ta1,tb1,tup,tdup,tun,tdun), ut(tp,ta2,tb2,tup,tdup,tun,tdun))
print('F1 scores:')
c(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
c(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
c(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
c(acc(tp,ta1,tb1), acc(tp,ta2,tb2))




tp <- 0.6
ta1 <- 0.6
tb1 <- 0.7
ta2 <- 0.9
tb2 <- 0.6
tup <- 15 -2
tdup <- tup + 5 +2
tun <- 25 -2
tdun <- tun + 75+2
##
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
matrix(c(ta1*tp,(1-ta1)*tp,
         (1-tb1)*(1-tp),tb1*(1-tp)), 2,2)
print('CM2:')
matrix(c(ta2*tp,(1-ta2)*tp,
         (1-tb2)*(1-tp),tb2*(1-tp)), 2,2)
print('UM:')
matrix(c(tup, tup-tdup,
         tun-tdun, tun),2,2)
##
print('utilities:')
c(ut(tp,ta1,tb1,tup,tdup,tun,tdun), ut(tp,ta2,tb2,tup,tdup,tun,tdun))
print('F1 scores:')
c(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
c(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
c(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
c(acc(tp,ta1,tb1), acc(tp,ta2,tb2))


tp <- 1-0.8
ta1 <- 0.9
tb1 <- 0.6
ta2 <- 0.8
tb2 <- 0.9
tun <- 0
tdun <- tun + 1
tup <- 15
tdup <- tup + 80
##
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
matrix(c(ta1*tp,(1-ta1)*tp,
         (1-tb1)*(1-tp),tb1*(1-tp)), 2,2)
print('CM2:')
matrix(c(ta2*tp,(1-ta2)*tp,
         (1-tb2)*(1-tp),tb2*(1-tp)), 2,2)
print('UM:')
matrix(c(tup, tup-tdup,
         tun-tdun, tun),2,2)
##
print('utilities:')
c(ut(tp,ta1,tb1,tup,tdup,tun,tdun), ut(tp,ta2,tb2,tup,tdup,tun,tdun))
print('F1 scores:')
c(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
c(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))

############################################################

set.seed(149)
nn <- 10^6
lp <- runif(nn, 0.5, 1)
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
lup <- runif(nn, 0, 1)
ldup <- runif(nn, 0, 1)
lun <- runif(nn, 0, 1)
ldun <- runif(nn, 0, 1)
##
ut1 <- ut(lp,la1,lb1,lup,ldup,lun,ldun)
ut2 <- ut(lp,la2,lb2,lup,ldup,lun,ldun)
Du <- ut2-ut1
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)

tplot(x=Du, y=Df,
      type='p',pch='.', xlab='U(1)-U(2)', ylab='F1(1)-F1(2)')

tplot(x=Du, y=Dm,
      type='p',pch='.', xlab='U(1)-U(2)', ylab='MCC(1)-MCC(2)')

sele <- which( (Du > 0 & Df < 0 & Dm < 0) | (Du < 0 & Df > 0 & Dm > 0) )

sele <- which( (ut2 > 0 & ut1 < 0 & Df < 0 & Dm < 0) | (ut2 < 0 & ut1 > 0 & Df > 0 & Dm > 0) )

tplot(x=Du[sele], y=Df[sele],
      type='p',pch='.', xlab='U(1)-U(2)', ylab='F1(1)-F1(2)')

tplot(x=Du[sele], y=Dm[sele],
      type='p',pch='.', xlab='U(1)-U(2)', ylab='MCC(1)-MCC(2)')

## okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele]) +ldup[sele]+ldun[sele])]
okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele]) +ldup[sele]+ldun[sele]+lp[sele])]
##
print('P+:')
lp[okpoint]
c(la1[okpoint], lb1[okpoint])
c(la2[okpoint], lb2[okpoint])
##
print('CM1:')
matrix(c(la1[okpoint]*lp[okpoint],(1-la1[okpoint])*lp[okpoint],
         (1-lb1[okpoint])*(1-lp[okpoint]),lb1[okpoint]*(1-lp[okpoint])), 2,2)
print('CM2:')
matrix(c(la2[okpoint]*lp[okpoint],(1-la2[okpoint])*lp[okpoint],
         (1-lb2[okpoint])*(1-lp[okpoint]),lb2[okpoint]*(1-lp[okpoint])), 2,2)
print('UM:')
matrix(c(lup[okpoint], lup[okpoint]-ldup[okpoint],
         lun[okpoint]-ldun[okpoint], lun[okpoint]),2,2)
##
print('utilities:')
c(ut(lp,la1,lb1,lup,ldup,lun,ldun)[okpoint], ut(lp,la2,lb2,lup,ldup,lun,ldun)[okpoint])
print('F1 scores:')
c(f1score(lp,la1,lb1)[okpoint], f1score(lp,la2,lb2)[okpoint])
print('MCCs:')
c(mcc(lp,la1,lb1)[okpoint], mcc(lp,la2,lb2)[okpoint])


tp <- 0.8
ta1 <- 0.6
tb1 <- 0.9
ta2 <- 0.9
tb2 <- 0.8
tup <- 0
tdup <- tup + 1
tun <- 15
tdun <- tun + 80
##
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
matrix(c(ta1*tp,(1-ta1)*tp,
         (1-tb1)*(1-tp),tb1*(1-tp)), 2,2)
print('CM2:')
matrix(c(ta2*tp,(1-ta2)*tp,
         (1-tb2)*(1-tp),tb2*(1-tp)), 2,2)
print('UM:')
matrix(c(tup, tup-tdup,
         tun-tdun, tun),2,2)
##
print('utilities:')
c(ut(tp,ta1,tb1,tup,tdup,tun,tdun), ut(tp,ta2,tb2,tup,tdup,tun,tdun))
print('F1 scores:')
c(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
c(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
c(ta1*tp/(ta1*tp+(1-tb1)*(1-tp)), ta2*tp/(ta2*tp+(2-tb1)*(1-tp)))



######################################################

set.seed(149)
nn <- 10^7
lp <- runif(nn, 0.5, 1)
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
lutp <- runif(nn, 0, 1)
lufp <- runif(nn, -1, 0)
lutn <- runif(nn, 0, 1)
lufn <- runif(nn, -1, 0)
##
ut1 <- ut(lp,la1,lb1,lutp,lutn,lufp,lufn)
ut2 <- ut(lp,la2,lb2,lutp,lutn,lufp,lufn)
Du <- ut2-ut1
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
Dac <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)

tplot(x=Du, y=Df,
      type='p',pch='.', xlab='U(1)-U(2)', ylab='F1(1)-F1(2)')

tplot(x=Du, y=Dm,
      type='p',pch='.', xlab='U(1)-U(2)', ylab='MCC(1)-MCC(2)')


## sele <- which( (ut2 > 0 & ut1 < 0 & Df < 0 & Dm < 0) | (ut2 < 0 & ut1 > 0 & Df > 0 & Dm > 0) )

set.seed(149)
nn <- 10^6
lp <- runif(nn, 0.5, 1)
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
lutp <- runif(nn, 0, 1)
lufp <- runif(nn, -1, 0)
lutn <- runif(nn, 0, 1)
lufn <- runif(nn, -1, 0)
##
ut1 <- ut(lp,la1,lb1,lutp,lutn,lufp,lufn)
ut2 <- ut(lp,la2,lb2,lutp,lutn,lufp,lufn)
Du <- ut2-ut1
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
Dac <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)
## sele <- which( (ut2 > 0 & ut1 < 0 & Df < 0 & Dm < 0 & Dpr < 0 & Dac < 0) | (ut2 < 0 & ut1 > 0 & Df > 0 & Dm > 0 & Dpr > 0 & Dac > 0) )
sele <- which(Du < 0 & Df > 0 & Dm > 0 & Dpr > 0 & Dac > 0 & lutn > ut2)# & ut2 >= lutp)
length(sele)
## okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele]) +ldup[sele]+ldun[sele])]
## okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele])+abs(Dpr[sele])+abs(Dac[sele]) +ldup[sele]+ldun[sele]+lp[sele])]
okpoint <- sele[which.max(abs(Du[sele])*3+abs(Df[sele])+abs(Dm[sele])+abs(Dpr[sele])+abs(Dac[sele]) -2*lp[sele])]
##
lp <- lp[okpoint]
la1 <- la1[okpoint]
lb1 <- lb1[okpoint]
la2 <- la2[okpoint]
lb2 <- lb2[okpoint]
lutp <- lutp[okpoint]
lutn <- lutn[okpoint]
lufp <- lufp[okpoint]
lufn <- lufn[okpoint]
##

Su <- -1
Xu <- 100
##
print('---------------------------------')
print('P+:')
lp
c(la1, lb1)
c(la2, lb2)
##
print('CM1:')
matrix(c(la1*lp,(1-la1)*lp,
         (1-lb1)*(1-lp),lb1*(1-lp)), 2,2)
print('CM2:')
matrix(c(la2*lp,(1-la2)*lp,
         (1-lb2)*(1-lp),lb2*(1-lp)), 2,2)
print('UM:')
Su+Xu*matrix(c(lutp, lufp,
         lufn, lutn),2,2)
##
print('utilities:')
Su+Xu*c(ut(lp,la1,lb1,lutp,lutn,lufp,lufn), ut(lp,la2,lb2,lutp,lutn,lufp,lufn))
print('F1 scores:')
c(f1score(lp,la1,lb1), f1score(lp,la2,lb2))
print('MCCs:')
c(mcc(lp,la1,lb1), mcc(lp,la2,lb2))
print('Precs:')
c(prec(lp,la1,lb1), prec(lp,la2,lb2))
print('Accs:')
c(acc(lp,la1,lb1), acc(lp,la2,lb2))

##########################################################
#### candidate for example in paper
##########################################################
tp <- 0.5
ta1 <- 0.74
tb1 <- 0.28*2
ta2 <- 0.60
tb2 <- 0.49*2
tutp <- 60
tufp <- 0
tufn <- -85
tutn <- 5
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
matrix(c(ta1*tp,(1-ta1)*tp, (1-tb1)*(1-tp),tb1*(1-tp)), 2,2)
print('CM2:')
matrix(c(ta2*tp,(1-ta2)*tp, (1-tb2)*(1-tp),tb2*(1-tp)), 2,2)
print('UM:')
matrix(c(tutp, tufn, tufp, tutn),2,2)
##
print('utilities:')
c(ut(tp,ta1,tb1,tutp,tutn,tufp,tufn), ut(tp,ta2,tb2,tutp,tutn,tufp,tufn))
print('F1 scores:')
c(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
c(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
c(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
c(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('AUCs:')
c(auc(tp,ta1,tb1), auc(tp,ta2,tb2))
##########################################################
##########################################################





##########################################################
#### With fixed utility matrix
##########################################################
set.seed(149)
nn <- 10^4
lp <- runif(nn, 0.5, 1)
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
lutp <- 1
lutn <- 1
lufp <- -10
lufn <- 0
##
ut1 <- ut(lp,la1,lb1,lutp,lutn,lufp,lufn)
ut2 <- ut(lp,la2,lb2,lutp,lutn,lufp,lufn)
Du <- ut2-ut1
f1 <- f1score(lp,la1,lb1)
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
m1 <- mcc(lp,la1,lb1)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
pr1 <- prec(lp,la1,lb1)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
ac1 <- acc(lp,la1,lb1)
Dac <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)
bac1 <- bacc(lp,la1,lb1)
pdff(paste0('example_discrepancies_',lutp,'_',lutn,'_',lufp,'_',lufn))
par(mfcol=c(2,3))
tplot(x=c(0,1),y=c(0,1),col='white',xticks=NA,yticks=NA,xlab=NA, ylab=NA)
text(x=0.5,y=0.5,labels=paste0('$TP = ', lutp, '  $FP = ', lufp,
                               '\n$FN = ', lufn, '  $TN = ', lutn),
     cex=2.5, family='Palatino')
#tplot(x=c(0,1),y=c(0,1),col='white',xticks=NA,yticks=NA,xlab=NA, ylab=NA)
tplot(x=ut1, y=f1, type='p', pch='.', xlab='utility', ylab='F1 score', family='Palatino')
tplot(x=ut1, y=m1, type='p', pch='.', xlab='utility', ylab='MCC', family='Palatino')
tplot(x=ut1, y=pr1, type='p', pch='.', xlab='utility', ylab='Precision', family='Palatino')
tplot(x=ut1, y=ac1, type='p', pch='.', xlab='utility', ylab='Accuracy', family='Palatino')
tplot(x=ut1, y=bac1, type='p', pch='.', xlab='utility', ylab='Bal. accuracy', family='Palatino')
dev.off()



## sele <- which( (ut2 > 0 & ut1 < 0 & Df < 0 & Dm < 0 & Dpr < 0 & Dac < 0) | (ut2 < 0 & ut1 > 0 & Df > 0 & Dm > 0 & Dpr > 0 & Dac > 0) )
sele <- which(Du<0 & Dpr>0 & Df>0 & Dm>0 & Dac>0)
length(sele)
## okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele]) +ldup[sele]+ldun[sele])]
## okpoint <- sele[which.max(abs(Du[sele])+abs(Df[sele])+abs(Dm[sele])+abs(Dpr[sele])+abs(Dac[sele]) +ldup[sele]+ldun[sele]+lp[sele])]
okpoint <- sele[which.max(abs(Du[sele])*100+abs(Df[sele])+abs(Dm[sele])+abs(Dpr[sele])+abs(Dac[sele]) +la1[sele]+lb1[sele]+la2[sele]+lb2[sele])]
##
lp <- lp[okpoint]
la1 <- la1[okpoint]
lb1 <- lb1[okpoint]
la2 <- la2[okpoint]
lb2 <- lb2[okpoint]
##
Su <- 0
Xu <- 1
##
print('---------------------------------')
print('P+:')
lp
c(la1, lb1)
c(la2, lb2)
##
print('CM1:')
matrix(c(la1*lp,(1-la1)*lp, (1-lb1)*(1-lp),lb1*(1-lp)), 2,2)
print('CM2:')
matrix(c(la2*lp,(1-la2)*lp, (1-lb2)*(1-lp),lb2*(1-lp)), 2,2)
print('UM:')
Su+Xu*matrix(c(lutp, lufn, lufp, lutn),2,2)
##
print('utilities:')
Su+Xu*c(ut(lp,la1,lb1,lutp,lutn,lufp,lufn), ut(lp,la2,lb2,lutp,lutn,lufp,lufn))
print('F1 scores:')
c(f1score(lp,la1,lb1), f1score(lp,la2,lb2))
print('MCCs:')
c(mcc(lp,la1,lb1), mcc(lp,la2,lb2))
print('Precs:')
c(prec(lp,la1,lb1), prec(lp,la2,lb2))
print('Accs:')
c(acc(lp,la1,lb1), acc(lp,la2,lb2))




##########################################################
#### With fixed utility matrix & P-probability
##########################################################
set.seed(149)
nn <- 10^4
pdff(paste0('example_discrepancies_P_',lutp,'_',lutn,'_',lufp,'_',lufn))
for(i in c(0.5, 0.75, 0.9, 0.99, 0.999)){
lp <- i # runif(nn, 0.5, 1)
##
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
lutp <- 10
lutn <- 1
lufp <- 0
lufn <- 0
##
ut1 <- ut(lp,la1,lb1,lutp,lutn,lufp,lufn)
ut2 <- ut(lp,la2,lb2,lutp,lutn,lufp,lufn)
Du <- ut2-ut1
f1 <- f1score(lp,la1,lb1)
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
m1 <- mcc(lp,la1,lb1)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
pr1 <- prec(lp,la1,lb1)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
ac1 <- acc(lp,la1,lb1)
Dac <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)
bac1 <- bacc(lp,la1,lb1)
par(mfcol=c(2,3))
tplot(x=c(0,1),y=c(0,1),col='white',xticks=NA,yticks=NA,xlab=NA, ylab=NA)
text(x=0.5,y=0.5,labels=paste0('$TP = ', lutp, '  $FP = ', lufp,
                               '\n$FN = ', lufn, '  $TN = ', lutn,
                               '\nPr(P) = ',lp),
     cex=2.5, family='Palatino')
#tplot(x=c(0,1),y=c(0,1),col='white',xticks=NA,yticks=NA,xlab=NA, ylab=NA)
tplot(x=ut1, y=f1, type='p', pch='.', xlab='utility', ylab='F1 score', family='Palatino')
tplot(x=ut1, y=m1, type='p', pch='.', xlab='utility', ylab='MCC', family='Palatino')
tplot(x=ut1, y=pr1, type='p', pch='.', xlab='utility', ylab='Precision', family='Palatino')
tplot(x=ut1, y=ac1, type='p', pch='.', xlab='utility', ylab='Accuracy', family='Palatino')
tplot(x=ut1, y=bac1, type='p', pch='.', xlab='utility', ylab='Bal. accuracy', family='Palatino')
}
dev.off()

###########################################################
#### With various utility matrices and positive-class probs
###########################################################
set.seed(149)
nn <- 10^5
lp <- runif(nn, 0.5, 1)
##
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
basetp <- matrix(c(1,0,0,0),2,2)
basetn <- matrix(c(0,0,0,1),2,2)
basefpfn <- array(c(1/2,1/2,0,0,  0,0,1/2,1/2),dim=c(2,2,2))
sides <- sample(1:2,nn,replace=T)
convpoints <- LaplacesDemon::rdirichlet(n=nn, alpha=rep(1,3))
lum <- sapply(1:nn, function(i){
    temp <- basetp * convpoints[i,1] + basetn * convpoints[i,2] +
        basefpfn[,,sides[i]] * convpoints[i,3]
    temp <- temp - min(temp)
    temp <- temp/max(temp)
})
dim(lum) <- c(2,2,nn)
## testlum <- lapply(1:nn, function(i){
##     basetp * convpoints[i,1] + basetn * convpoints[i,2] +
##         basefpfn[,,sides[i]] * convpoints[i,3]
## }
## )
##
percerror <- 1/4
lumappr <- sapply(1:nn, function(i){
    temp <- lum[,,i] + matrix(rnorm(4),2,2)*lum[,,i]*percerror
})
dim(lumappr) <- c(2,2,nn)
##
ut1 <- sapply(1:nn, function(i){ut(lp[i], la1[i], lb1[i], lum[,,i])})
ut2 <- sapply(1:nn, function(i){ut(lp[i], la2[i], lb2[i], lum[,,i])})
Du <- ut2-ut1
#f1 <- f1score(lp,la1,lb1)
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
#m1 <- mcc(lp,la1,lb1)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
#pr1 <- prec(lp,la1,lb1)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
#ac1 <- acc(lp,la1,lb1)
Dacc <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)
Dbacc <- -bacc(lp,la1,lb1)+bacc(lp,la2,lb2)
##
Duappr <- sapply(1:nn, function(i){
    ut(lp[i], la2[i], lb2[i], lumappr[,,i]) - ut(lp[i], la1[i], lb1[i], lumappr[,,i])
})
##
pdff(paste0('example_discrepancies_all'))
posi <- (Du > 0 & Df > 0) | (Du < 0 & Df < 0)
tplot(x=Du[posi], y=Df[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. F1 score', family='Palatino')
tplot(x=Du[!posi], y=Df[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. F1 score', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Dm > 0) | (Du < 0 & Dm < 0)
tplot(x=Du[posi], y=Dm[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. MCC', family='Palatino')
tplot(x=Du[!posi], y=Dm[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. MCC', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Dpr > 0) | (Du < 0 & Dpr < 0)
tplot(x=Du[posi], y=Dpr[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. Precision', family='Palatino')
tplot(x=Du[!posi], y=Dpr[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. Precision', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Dbacc > 0) | (Du < 0 & Dbacc < 0)
tplot(x=Du[posi], y=Dbacc[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. Bal. accuracy', family='Palatino')
tplot(x=Du[!posi], y=Dbacc[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. Bal. accuracy', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Dacc > 0) | (Du < 0 & Dacc < 0)
tplot(x=Du[posi], y=Dacc[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. Accuracy', family='Palatino')
tplot(x=Du[!posi], y=Dacc[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. Accuracy', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Duappr > 0) | (Du < 0 & Duappr < 0)
tplot(x=Du[posi], y=Duappr[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab=paste0('diff. utility with ',round(100*percerror),'% error'), family='Palatino')
tplot(x=Du[!posi], y=Duappr[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab=paste0('diff. ',round(100*percerror),'% approx utility'), family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
dev.off()





###########################################################
#### With various utility matrices and positive-class probs
#### case for "better" confusion matrices
###########################################################
set.seed(149)
nn <- 10^5
lp <- runif(nn, 0.5, 1)
##
la1 <- 0.5+0.5*rbeta(10^5, shape1=3, shape2=1)
lb1 <- 0.5+0.5*rbeta(10^5, shape1=3, shape2=1)
la2 <- 0.5+0.5*rbeta(10^5, shape1=3, shape2=1)
lb2 <- 0.5+0.5*rbeta(10^5, shape1=3, shape2=1)
basetp <- matrix(c(1,0,0,0),2,2)
basetn <- matrix(c(0,0,0,1),2,2)
basefpfn <- array(c(1/2,1/2,0,0,  0,0,1/2,1/2),dim=c(2,2,2))
sides <- sample(1:2,nn,replace=T)
convpoints <- LaplacesDemon::rdirichlet(n=nn, alpha=rep(1,3))
lum <- sapply(1:nn, function(i){
    temp <- basetp * convpoints[i,1] + basetn * convpoints[i,2] +
        basefpfn[,,sides[i]] * convpoints[i,3]
    temp <- temp - min(temp)
    temp <- temp/max(temp)
})
dim(lum) <- c(2,2,nn)
## testlum <- lapply(1:nn, function(i){
##     basetp * convpoints[i,1] + basetn * convpoints[i,2] +
##         basefpfn[,,sides[i]] * convpoints[i,3]
## }
## )
##
percerror <- 1/4
lumappr <- sapply(1:nn, function(i){
    temp <- lum[,,i] + matrix(rnorm(4),2,2)*lum[,,i]*percerror
})
dim(lumappr) <- c(2,2,nn)
##
ut1 <- sapply(1:nn, function(i){ut(lp[i], la1[i], lb1[i], lum[,,i])})
ut2 <- sapply(1:nn, function(i){ut(lp[i], la2[i], lb2[i], lum[,,i])})
Du <- ut2-ut1
#f1 <- f1score(lp,la1,lb1)
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
#m1 <- mcc(lp,la1,lb1)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
#pr1 <- prec(lp,la1,lb1)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
#ac1 <- acc(lp,la1,lb1)
Dacc <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)
Dbacc <- -bacc(lp,la1,lb1)+bacc(lp,la2,lb2)
##
Duappr <- sapply(1:nn, function(i){
    ut(lp[i], la2[i], lb2[i], lumappr[,,i]) - ut(lp[i], la1[i], lb1[i], lumappr[,,i])
})
##
pdff(paste0('example_discrepancies_better'))
posi <- (Du > 0 & Df > 0) | (Du < 0 & Df < 0)
tplot(x=Du[posi], y=Df[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. F1 score', family='Palatino')
tplot(x=Du[!posi], y=Df[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. F1 score', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Dm > 0) | (Du < 0 & Dm < 0)
tplot(x=Du[posi], y=Dm[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. MCC', family='Palatino')
tplot(x=Du[!posi], y=Dm[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. MCC', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Dpr > 0) | (Du < 0 & Dpr < 0)
tplot(x=Du[posi], y=Dpr[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. Precision', family='Palatino')
tplot(x=Du[!posi], y=Dpr[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. Precision', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Dbacc > 0) | (Du < 0 & Dbacc < 0)
tplot(x=Du[posi], y=Dbacc[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. Bal. accuracy', family='Palatino')
tplot(x=Du[!posi], y=Dbacc[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. Bal. accuracy', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Dacc > 0) | (Du < 0 & Dacc < 0)
tplot(x=Du[posi], y=Dacc[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab='diff. Accuracy', family='Palatino')
tplot(x=Du[!posi], y=Dacc[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab='diff. Accuracy', family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
posi <- (Du > 0 & Duappr > 0) | (Du < 0 & Duappr < 0)
tplot(x=Du[posi], y=Duappr[posi], type='p', pch='.', col=1, xlab='diff. utility', ylab=paste0('diff. utility with ',round(100*percerror),'% error'), family='Palatino')
tplot(x=Du[!posi], y=Duappr[!posi], type='p', pch='.', col=2, xlab='diff. utility', ylab=paste0('diff. ',round(100*percerror),'% approx utility'), family='Palatino', add=T)
legend('topleft',legend=paste0('correctly ranked cases: ',round(sum(posi)/nn*100),'%'), bty='n', cex=1.5)
##
dev.off()


###########################################################
#### Search for matrices for paper's example 2
###########################################################
set.seed(941)
nn <- 10^7
lp <- runif(nn, 0, 1)
##
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
lcm1 <- confm(lp,la1,lb1)
lcm2 <- confm(lp,la2,lb2)
##
basetpn <- matrix(c(1,0,0,1),2,2)
basefpfn <- array(c(1,1,0,1,  1,0,1,1),dim=c(2,2,2))
sides <- sample(1:2,nn,replace=T)
convpoints <- LaplacesDemon::rdirichlet(n=nn, alpha=rep(1,2))
## system.time(lum1 <- sapply(1:nn,function(i){
## #matrix(c(10,1,-10,1),2,2)
##     temp <- basetp * convpoints[i,1] + basetn * convpoints[i,2] +
##         basefpfn[,,sides[i]] * convpoints[i,3]
##     ## temp <- temp - min(temp)
##     temp <- temp/max(temp)
## }))
## dim(lum1) <- c(2,2,nn)
lum <- array(rep(c(10,0,11,12), nn), dim=c(2,2,nn))
## lum <- future_sapply(1:nn,function(i){
##     temp <- basetpn * convpoints[i,1] +
##         basefpfn[,,sides[i]] * convpoints[i,2]
##     ## ## temp <- temp - min(temp)
##     temp <- temp/max(temp)
## })
## dim(lum) <- c(2,2,nn)
## testlum <- lapply(1:nn, function(i){
##     basetp * convpoints[i,1] + basetn * convpoints[i,2] +
##         basefpfn[,,sides[i]] * convpoints[i,3]
## }
## )
##
## percerror <- 1/4
## lumappr <- sapply(1:nn, function(i){
##     temp <- lum[,,i] + matrix(rnorm(4),2,2)*lum[,,i]*percerror
## })
## dim(lumappr) <- c(2,2,nn)
##
ut1 <- ut(lcm1, lum)
ut2 <- ut(lcm2, lum)
Du <- ut2-ut1
#f1 <- f1score(lp,la1,lb1)
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
#m1 <- mcc(lp,la1,lb1)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
#pr1 <- prec(lp,la1,lb1)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
#ac1 <- acc(lp,la1,lb1)
Dacc <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)
Dbacc <- -bacc(lp,la1,lb1)+bacc(lp,la2,lb2)
##
## Duappr <- sapply(1:nn, function(i){
##     ut(lp[i], la2[i], lb2[i], lumappr[,,i]) - ut(lp[i], la1[i], lb1[i], lumappr[,,i])
## })
##
Drec <- la2-la1
Dspec <- lb2-lb1
## ##
## dcosts <- apply(lum,3,function(x){min(x[1,1],x[2,2])-max(x[1,2],x[2,1])})
## ##
minp <- future_sapply(1:nn,function(i){
    c(min(lum[1,1,i],lum[2,2,i]), max(lum[1,2,i],lum[2,1,i]))
})
##
ddiff <- minp[1,]-minp[2,]
##
dshifts <- ut1 - minp[2,]
## ##
dshifts2 <- minp[1,] - ut2
##
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0 & Drec<=0 & Dspec<=0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Drec>=0 & Dspec>=0)
## sum(select)
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0)
## sum(select)
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0 & Drec<=0 & Dspec<=0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Drec>=0 & Dspec>=0)
## sum(select)
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0 & Drec<0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Drec>0)
## sum(select)
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0 & Dspec<=0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Dspec>=0)
## sum(select)
##
select1 <- (2*Du/(ut1+ut2) < -0.1 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Drec>0)# & ddiff>0)
sum(select1)
select2 <- (2*Du/(ut1+ut2) < -0.1 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Dspec>0)# & ddiff>0)
sum(select2)

##
okpoint <- which.max(abs(Du[select])*0+abs(Df[select])+abs(Dm[select])+abs(Dpr[select])+abs(Dacc[select])+abs(Drec[select])-10*2*((Drec[select]>0)*abs(-lb1[select]+lb2[select])/(lb1[select]+lb2[select])+(Dspec[select]>0)*abs(-la1[select]+la2[select])/(la1[select]+la2[select]))-0*lp[select])

##
okpoint <- which.max(abs(Du[select])*0+abs(Df[select])+abs(Dm[select])+abs(Dpr[select])+abs(Dacc[select])+abs(Drec[select])-50*2*((Drec[select]>0)*abs(-lb1[select]+lb2[select])+(Dspec[select]>0)*abs(-la1[select]+la2[select]))-0*lp[select])

okpoint <- which.max(abs(Du[select])*0+abs(Df[select])+abs(Dm[select])+abs(Dpr[select])+abs(Dacc[select])+abs(Drec[select])+1000*2*((Drec[select]>0)*abs(lb2[select])+(Dspec[select]>0)*abs(la2[select]))-0*lp[select])

##
select <- select2
okpoint <- which.max(abs(Du[select])*0+abs(Df[select])+abs(Dm[select])+abs(Dpr[select])+abs(Dacc[select])+abs(Drec[select])-10*abs(-la1[select]+la2[select])/(la1[select]+la2[select])-0*lp[select])
#okpoint <- 972#615#292#145
sigf <- 3
tp <- signif(lp[select][okpoint] ,sigf)
ta1 <- la1[select][okpoint]
tb1 <- lb1[select][okpoint]
tcm1 <- round(confm(tp,ta1,tb1),sigf)
ta2 <- la2[select][okpoint]
tb2 <- lb2[select][okpoint]
tcm2 <- round(confm(tp,ta2,tb2),sigf)
##
Xu <- 1
Su <- 0# round(-(tut1*1/4+tut2*3/4)*Xu, 0)
tum <- signif(round(Su+Xu*lum[,,select][,,okpoint] ,0),sigf)
#tum <- matrix(c(10,0,-600,450),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),3)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)


###########################################################
#### Search for matrices for paper's example
###########################################################
set.seed(941)
nn <- 10^7
lp <- runif(nn, 0.5, 1)
##
la1 <- runif(nn, 0.5, 1)
lb1 <- runif(nn, 0.5, 1)
la2 <- runif(nn, 0.5, 1)
lb2 <- runif(nn, 0.5, 1)
lcm1 <- confm(lp,la1,lb1)
lcm2 <- confm(lp,la2,lb2)
##
basetp <- matrix(c(1,0,0,0),2,2)
basetn <- matrix(c(0,0,0,1),2,2)
basefpfn <- array(c(1/2,1/2,0,0,  0,0,1/2,1/2),dim=c(2,2,2))
sides <- sample(1:2,nn,replace=T)
convpoints <- LaplacesDemon::rdirichlet(n=nn, alpha=rep(1,3))
## system.time(lum1 <- sapply(1:nn,function(i){
## #matrix(c(10,1,-10,1),2,2)
##     temp <- basetp * convpoints[i,1] + basetn * convpoints[i,2] +
##         basefpfn[,,sides[i]] * convpoints[i,3]
##     ## temp <- temp - min(temp)
##     temp <- temp/max(temp)
## }))
## dim(lum1) <- c(2,2,nn)
system.time(lum <- future_sapply(1:nn,function(i){
#matrix(c(10,1,-10,1),2,2)
    temp <- basetp * convpoints[i,1] + basetn * convpoints[i,2] +
        basefpfn[,,sides[i]] * convpoints[i,3]
    ## temp <- temp - min(temp)
    temp <- temp/max(temp)
}))
dim(lum) <- c(2,2,nn)
## testlum <- lapply(1:nn, function(i){
##     basetp * convpoints[i,1] + basetn * convpoints[i,2] +
##         basefpfn[,,sides[i]] * convpoints[i,3]
## }
## )
##
## percerror <- 1/4
## lumappr <- sapply(1:nn, function(i){
##     temp <- lum[,,i] + matrix(rnorm(4),2,2)*lum[,,i]*percerror
## })
## dim(lumappr) <- c(2,2,nn)
##
ut1 <- ut(lcm1, lum)
ut2 <- ut(lcm2, lum)
Du <- ut2-ut1
#f1 <- f1score(lp,la1,lb1)
Df <- -f1score(lp,la1,lb1)+f1score(lp,la2,lb2)
#m1 <- mcc(lp,la1,lb1)
Dm <- -mcc(lp,la1,lb1)+mcc(lp,la2,lb2)
#pr1 <- prec(lp,la1,lb1)
Dpr <- -prec(lp,la1,lb1)+prec(lp,la2,lb2)
#ac1 <- acc(lp,la1,lb1)
Dacc <- -acc(lp,la1,lb1)+acc(lp,la2,lb2)
Dbacc <- -bacc(lp,la1,lb1)+bacc(lp,la2,lb2)
##
## Duappr <- sapply(1:nn, function(i){
##     ut(lp[i], la2[i], lb2[i], lumappr[,,i]) - ut(lp[i], la1[i], lb1[i], lumappr[,,i])
## })
##
Drec <- la2-la1
Dspec <- lb2-lb1
## ##
## dcosts <- apply(lum,3,function(x){min(x[1,1],x[2,2])-max(x[1,2],x[2,1])})
## ##
minp <- future_sapply(1:nn,function(i){
    c(min(lum[1,1,i],lum[2,2,i]), max(lum[1,2,i],lum[2,1,i]))
})
##
ddiff <- minp[1,]-minp[2,]
##
dshifts <- ut1 - minp[2,]
## ##
dshifts2 <- minp[1,] - ut2


## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0 & Drec<=0 & Dspec<=0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Drec>=0 & Dspec>=0)
## sum(select)
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0)
## sum(select)
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0 & Drec<=0 & Dspec<=0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Drec>=0 & Dspec>=0)
## sum(select)
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0 & Drec<0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Drec>0)
## sum(select)
## select <- (Du>0 & Df<0 & Dm<0 & Dpr<0 & Dacc<0 & Dbacc<0 & Dspec<=0) |
##     (Du<0 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & Dspec>=0)
## sum(select)


select <- (2*Du/(ut1+ut2) < -0.1 & Df>0 & Dm>0 & Dpr>0 & Dacc>0 & Dbacc>0 & (Drec>0 | Dspec>0))# & ddiff>0)
sum(select)

##
okpoint <- which.max(abs(Du[select])*0+abs(Df[select])+abs(Dm[select])+abs(Dpr[select])+abs(Dacc[select])+abs(Drec[select])-10*2*((Drec[select]>0)*abs(-lb1[select]+lb2[select])/(lb1[select]+lb2[select])+(Dspec[select]>0)*abs(-la1[select]+la2[select])/(la1[select]+la2[select]))-0*lp[select])

##
okpoint <- which.max(abs(Du[select])*0+abs(Df[select])+abs(Dm[select])+abs(Dpr[select])+abs(Dacc[select])+abs(Drec[select])-10*2*((Drec[select]>0)*abs(-lb1[select]+lb2[select])+(Dspec[select]>0)*abs(-la1[select]+la2[select]))-0*lp[select])

okpoint <- which.max(abs(Du[select])*0+abs(Df[select])+abs(Dm[select])+abs(Dpr[select])+abs(Dacc[select])+abs(Drec[select])+1000*2*((Drec[select]>0)*abs(lb2[select])+(Dspec[select]>0)*abs(la2[select]))-0*lp[select])

##okpoint <- 2649
sigf <- 3
tp <- signif(lp[select][okpoint] ,sigf)
ta1 <- la1[select][okpoint]
tb1 <- lb1[select][okpoint]
tcm1 <- round(confm(tp,ta1,tb1),sigf)
ta2 <- la2[select][okpoint]
tb2 <- lb2[select][okpoint]
tcm2 <- round(confm(tp,ta2,tb2),sigf)
##
Xu <- 1000
Su <- -139# round(-(tut1*1/4+tut2*3/4)*Xu, 0)
tum <- signif(round(Su+Xu*lum[,,select][,,okpoint] ,0),sigf)
#tum <- matrix(c(10,0,-600,450),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),3)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)

## [1] "P+:"
## > [1] 0.51
## > [1] 0.55 0.84
## > [1] 0.97 0.82
## > > [1] "CM1:"
## >        [,1]   [,2]
## [1,] 0.2805 0.0784
## [2,] 0.2295 0.4116
## > [1] "CM2:"
## >        [,1]   [,2]
## [1,] 0.4947 0.0882
## [2,] 0.0153 0.4018
## > [1] "UM:"
## >      [,1] [,2]
## [1,]    3  -59
## [2,]    0   41
## > > [1] "utilities:"
## > [1] 13.0915 12.7541
## > [1] "F1 scores:"
## > [1] 0.6456439 0.9052978
## > [1] "MCCs:"
## > [1] 0.4064416 0.8009273
## > [1] "Precs:"
## > [1] 0.7815548 0.8486876
## > [1] "Accs:"
## > [1] 0.6921 0.8965
## > [1] "Recalls:"
## > [1] 0.55 0.97
## > [1] "Specificities:"
## > [1] 0.84 0.82


tcm1 <- matrix(c(0.28,0.23,0.08,0.41),2,2)
tcm2 <- matrix(c(0.49,0.02,0.09,0.40),2,2)
tcoe1 <- cm2st(tcm1)
tcoe2 <- cm2st(tcm2)
tp <- tcoe1[1]
ta1 <- tcoe1[2]
tb1 <- tcoe1[3]
ta2 <- tcoe2[2]
tb2 <- tcoe2[3]
##
tum <- matrix(c(-105,-135,-625,275),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),4)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)


## > [1] "P+:"
## > [1] 0.88
## > [1] 0.5193395 0.8032112
## > [1] 0.9997283 0.7758837
## > > [1] "CM1:"
## >      [,1] [,2]
## [1,] 0.46 0.02
## [2,] 0.42 0.10
## > [1] "CM2:"
## >      [,1] [,2]
## [1,] 0.88 0.03
## [2,] 0.00 0.09
## > [1] "UM:"
## >      [,1] [,2]
## [1,]  -92  -95
## [2,]  -95  900
## > > > [1] "utilities:"
## > [1]    5.88   -2.81      NA -566.00
## > [1] "F1 scores:"
## > [1]  0.672  0.985     NA 37.800
## > [1] "MCCs:"
## > [1]   0.210   0.866      NA 122.000
## > [1] "Precs:"
## > [1] 0.951 0.970    NA 2.030
## > [1] "Accs:"
## > [1]  0.553  0.973     NA 55.000
## > [1] "Recalls:"
## > [1]  0.519  1.000     NA 63.200
## > [1] "Specificities:"
## > [1]  0.803  0.776     NA -3.460
## [1,]    3  -59
## [2,]    0   41
## > > [1] "utilities:"
## > [1] 13.0915 12.7541
## > [1] "F1 scores:"
## > [1] 0.6456439 0.9052978
## > [1] "MCCs:"
## > [1] 0.4064416 0.8009273
## > [1] "Precs:"
## > [1] 0.7815548 0.8486876
## > [1] "Accs:"
## > [1] 0.6921 0.8965
## > [1] "Recalls:"
## > [1] 0.55 0.97
## > [1] "Specificities:"
## > [1] 0.84 0.82


## > [1] "P+:"
## > [1] 0.7
## > [1] 0.5519250 0.5186422
## > [1] 0.5433391 0.9922608
## > > [1] "CM1:"
## >      [,1] [,2]
## [1,] 0.39 0.14
## [2,] 0.31 0.16
## > [1] "CM2:"
## >      [,1] [,2]
## [1,] 0.38  0.0
## [2,] 0.32  0.3
## > [1] "UM:"
## >      [,1] [,2]
## [1,]  340  240
## [2,] -660  260
## > > > [1] "utilities:"
## > [1]    3.2   -4.0     NA 1800.0
## > [1] "F1 scores:"
## > [1]  0.628  0.703     NA 11.200
## > [1] "MCCs:"
## > [1]   0.0648   0.5050       NA 155.0000
## > [1] "Precs:"
## > [1]  0.728  0.994     NA 30.900
## > [1] "Accs:"
## > [1]  0.542  0.678     NA 22.300
## > [1] "Recalls:"
## > [1]  0.552  0.543     NA -1.570
## > [1] "Specificities:"
## > [1]  0.519  0.992     NA 62.700


tcm1 <- matrix(c(0.39,0.31,0.14,0.16),2,2)
tcm2 <- matrix(c(0.38,0.32,0.0,0.3),2,2)
tcoe1 <- cm2st(tcm1)
tcoe2 <- cm2st(tcm2)
tp <- tcoe1[1]
ta1 <- tcoe1[2]
tb1 <- tcoe1[3]
ta2 <- tcoe2[2]
tb2 <- tcoe2[3]
##
tum <- matrix(c(132*1/4,-368*1/4,  202,222),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),4)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Bal. accs:')
com(bacc(tp,ta1,tb1), bacc(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)





tcm1 <- matrix(c(0.39,0.31,0.14,0.16),2,2)
tcm2 <- matrix(c(0.38,0.32,0.0,0.3),2,2)
tcoe1 <- cm2st(tcm1)
tcoe2 <- cm2st(tcm2)
tp <- tcoe1[1]
ta1 <- tcoe1[2]
tb1 <- tcoe1[3]
ta2 <- tcoe2[2]
tb2 <- tcoe2[3]
##
tum <- matrix(c(140,-375,200,220),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),4)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Bal. accs:')
com(bacc(tp,ta1,tb1), bacc(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)



tcm1 <- matrix(c(0.39,0.31,0.14,0.16),2,2)
tcm2 <- matrix(c(0.38,0.32,0.0,0.3),2,2)
tcoe1 <- cm2st(tcm1)
tcoe2 <- cm2st(tcm2)
tp <- tcoe1[1]
ta1 <- tcoe1[2]
tb1 <- tcoe1[3]
ta2 <- tcoe2[2]
tb2 <- tcoe2[3]
##
tum <- matrix(c(340,-660,240,260),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),4)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Bal. accs:')
com(bacc(tp,ta1,tb1), bacc(tp,ta2,tb2))
print('Krippendorff alphas:')
com(kri(tp,ta1,tb1), kri(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)


340,-660,240,260

tcm1 <- matrix(c(0.39,0.31,0.14,0.16),2,2)
tcm2 <- matrix(c(0.38,0.32,0.0,0.3),2,2)
tcoe1 <- cm2st(tcm1)
tcoe2 <- cm2st(tcm2)
tp <- tcoe1[1]
ta1 <- tcoe1[2]
tb1 <- tcoe1[3]
ta2 <- tcoe2[2]
tb2 <- tcoe2[3]
##
tum <- matrix(c(340,-660,240,260),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),4)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Bal. accs:')
com(bacc(tp,ta1,tb1), bacc(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)








tcm1 <- matrix(c(0.39,0.31,0.14,0.16),2,2)
tcm2 <- matrix(c(0.38,0.32,0.0,0.3),2,2)
tcoe1 <- cm2st(tcm1)
tcoe2 <- cm2st(tcm2)
tp <- tcoe1[1]
ta1 <- tcoe1[2]
tb1 <- tcoe1[3]
ta2 <- tcoe2[2]
tb2 <- tcoe2[3]
##
addt <- 30
tum <- matrix(c(340,-660,240-addt,260+addt),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),2)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Bal. accs:')
com(bacc(tp,ta1,tb1), bacc(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)





###############################################################
#### Example 2
###############################################################
## > [1] "P+:"
## > [1] 0.7
## > [1] 0.5571429 0.5333333
## > [1] 0.5428571 1.0000000
## > > [1] "CM1:"
## >      [,1] [,2]
## [1,] 0.39 0.14
## [2,] 0.31 0.16
## > [1] "CM2:"
## >      [,1] [,2]
## [1,] 0.38  0.0
## [2,] 0.32  0.3
## > [1] "UM:"
## >      [,1] [,2]
## [1,]  132  202
## [2,] -368  222
## > > > [1] "utilities:"
## > [1]     1.2    -1.0      NA -2200.0
## > [1] "F1 scores:"
## > [1]  0.6341  0.7037      NA 10.4000
## > [1] "MCCs:"
## > [1]   0.08307   0.51250        NA 144.20000
## > [1] "Precs:"
## > [1]  0.7358  1.0000      NA 30.4300
## > [1] "Accs:"
## > [1]  0.55  0.68    NA 21.14
## > [1] "Bal. accs:"
## > [1]  0.5452  0.7714      NA 34.3600
## > [1] "Recalls:"
## > [1]  0.5571  0.5429      NA -2.5970
## > [1] "Specificities:"
## > [1]  0.5333  1.0000      NA 60.8700
##
## > [1] "P+:"
## > [1] 0.7
## > [1] 0.5571429 0.5333333
## > [1] 0.5428571 1.0000000
## > > [1] "CM1:"
## >      [,1] [,2]
## [1,] 0.39 0.14
## [2,] 0.31 0.16
## > [1] "CM2:"
## >      [,1] [,2]
## [1,] 0.38  0.0
## [2,] 0.32  0.3
## > [1] "UM:"
## >      [,1] [,2]
## [1,]   33  202
## [2,]  -92  222
## > > > [1] "utilities:"
## > [1] 48.150 49.700     NA  3.168
## > [1] "F1 scores:"
## > [1]  0.6341  0.7037      NA 10.4000
## > [1] "MCCs:"
## > [1]   0.08307   0.51250        NA 144.20000
## > [1] "Precs:"
## > [1]  0.7358  1.0000      NA 30.4300
## > [1] "Accs:"
## > [1]  0.55  0.68    NA 21.14
## > [1] "Bal. accs:"
## > [1]  0.5452  0.7714      NA 34.3600
## > [1] "Recalls:"
## > [1]  0.5571  0.5429      NA -2.5970
## > [1] "Specificities:"
## > [1]  0.5333  1.0000      NA 60.8700

###############################################################
#### Example 1
###############################################################
## > [1] "P+:"
## > [1] 0.7
## > [1] 0.5571429 0.5333333
## > [1] 0.5428571 1.0000000
## > > [1] "CM1:"
## >      [,1] [,2]
## [1,] 0.39 0.14
## [2,] 0.31 0.16
## > [1] "CM2:"
## >      [,1] [,2]
## [1,] 0.38  0.0
## [2,] 0.32  0.3
## > [1] "UM:"
## >      [,1] [,2]
## [1,]  340  240
## [2,] -660  260
## > > > [1] "utilities:"
## > [1]    3.2   -4.0     NA 1800.0
## > [1] "F1 scores:"
## > [1]  0.6341  0.7037      NA 10.4000
## > [1] "MCCs:"
## > [1]   0.08307   0.51250        NA 144.20000
## > [1] "Precs:"
## > [1]  0.7358  1.0000      NA 30.4300
## > [1] "Accs:"
## > [1]  0.55  0.68    NA 21.14
## > [1] "Bal. accs:"
## > [1]  0.5452  0.7714      NA 34.3600
## > [1] "Krippendorff alphas:"
## > [1]  0.5249  0.6779      NA 25.4500
## > [1] "Recalls:"
## > [1]  0.5571  0.5429      NA -2.5970
## > [1] "Specificities:"
## > [1]  0.5333  1.0000      NA 60.8700
## small change:
## > [1] "P+:"
## > [1] 0.7
## > [1] 0.5571429 0.5333333
## > [1] 0.5428571 1.0000000
## > > [1] "CM1:"
## >      [,1] [,2]
## [1,] 0.39 0.14
## [2,] 0.31 0.16
## > [1] "CM2:"
## >      [,1] [,2]
## [1,] 0.38  0.0
## [2,] 0.32  0.3
## > [1] "UM:"
## >      [,1] [,2]
## [1,]  340  210
## [2,] -660  290
## > > > [1] "utilities:"
## > [1]  3.8  5.0   NA 27.0
## > [1] "F1 scores:"
## > [1]  0.63  0.70    NA 10.00
## > [1] "MCCs:"
## > [1]   0.083   0.510      NA 140.000
## > [1] "Precs:"
## > [1]  0.74  1.00    NA 30.00
## > [1] "Accs:"
## > [1]  0.55  0.68    NA 21.00
## > [1] "Bal. accs:"
## > [1]  0.55  0.77    NA 34.00
## > [1] "Recalls:"
## > [1]  0.56  0.54    NA -2.60
## > [1] "Specificities:"
## > [1]  0.53  1.00    NA 61.00





##
## > [1] "P+:"
## > [1] 0.7
## > [1] 0.5571429 0.5333333
## > [1] 0.5428571 1.0000000
## > > [1] "CM1:"
## >      [,1] [,2]
## [1,] 0.39 0.14
## [2,] 0.31 0.16
## > [1] "CM2:"
## >      [,1] [,2]
## [1,] 0.38  0.0
## [2,] 0.32  0.3
## > [1] "UM:"
## >      [,1] [,2]
## [1,]   85  240
## [2,] -165  260
## > > > [1] "utilities:"
## > [1] 57.2000 57.5000      NA  0.5231
## > [1] "F1 scores:"
## > [1]  0.6341  0.7037      NA 10.4000
## > [1] "MCCs:"
## > [1]   0.08307   0.51250        NA 144.20000
## > [1] "Precs:"
## > [1]  0.7358  1.0000      NA 30.4300
## > [1] "Accs:"
## > [1]  0.55  0.68    NA 21.14
## > [1] "Bal. accs:"
## > [1]  0.5452  0.7714      NA 34.3600
## > [1] "Recalls:"
## > [1]  0.5571  0.5429      NA -2.5970
## > [1] "Specificities:"
## > [1]  0.5333  1.0000      NA 60.8700



#########################################################
tcm1 <- matrix(c(0.39,0.31,0.14,0.16),2,2)
tcm2 <- matrix(c(0.38,0.32,0.0,0.3),2,2)
tcoe1 <- cm2st(tcm1)
tcoe2 <- cm2st(tcm2)
tp <- tcoe1[1]
ta1 <- tcoe1[2]
tb1 <- tcoe1[3]
ta2 <- tcoe2[2]
tb2 <- tcoe2[3]
##
set.seed(941)
nn <- 10^7
basetp <- matrix(c(1,0,0,0),2,2)
basetn <- matrix(c(0,0,0,1),2,2)
basefpfn <- array(c(1/2,1/2,0,0,  0,0,1/2,1/2),dim=c(2,2,2))
sides <- sample(1:2,nn,replace=T)
system.time(zum <- future_sapply(1:nn,function(i){
    temp <- basetp * convpoints[i,1] + basetn * convpoints[i,2] +
        basefpfn[,,sides[i]] * convpoints[i,3]
    temp <- temp/max(temp)
}))
dim(zum) <- c(2,2,nn)
upairs <- future_sapply(1:nn,function(i){
    c(ut(tcm1,zum[,,i]), ut(tcm2,zum[,,i]))
})


tplot(x=(upairs[1,1:10000]+upairs[2,1:10000])/2,y=upairs[2,1:10000]-upairs[1,1:10000],type='p',pch='.')


tum <- matrix(c(132*1/4,-368*1/4,  202,222),2,2)
tut1 <- ut(tcm1, tum)
tut2 <- ut(tcm2, tum)
##
print('---------------------------------')
print('P+:')
tp
c(ta1, tb1)
c(ta2, tb2)
##
print('CM1:')
tcm1
print('CM2:')
tcm2
print('UM:')
tum
##
com <- function(x,y){signif(c(x,y,NA,(y-x)/(x+y)*2*100),4)}
print('utilities:')
com(tut1, tut2)
print('F1 scores:')
com(f1score(tp,ta1,tb1), f1score(tp,ta2,tb2))
print('MCCs:')
com(mcc(tp,ta1,tb1), mcc(tp,ta2,tb2))
print('Precs:')
com(prec(tp,ta1,tb1), prec(tp,ta2,tb2))
print('Accs:')
com(acc(tp,ta1,tb1), acc(tp,ta2,tb2))
print('Bal. accs:')
com(bacc(tp,ta1,tb1), bacc(tp,ta2,tb2))
print('Recalls:')
com(ta1, ta2)
print('Specificities:')
com(tb1, tb2)



















