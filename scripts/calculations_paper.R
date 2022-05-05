## Author: PGL  Porta Mana
## Created: 2022-05-01T09:38:48+0200
## Last-Updated: 2022-05-04T20:00:32+0200
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
ut <- function(p,a,b,utm){
    sum(utm * matrix(c(p*a, p*(1-a), (1-p)*(1-b), (1-p)*b), 2,2))
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

