## Author: PGL  Porta Mana
## Created: 2022-04-15T11:48:45+0200
## Last-Updated: 2022-04-15T18:30:25+0200
################
## Calculation of joint probability for class & classifier-output
################
if(file.exists("/cluster/home/pglpm/R")){
    .libPaths(c("/cluster/home/pglpm/R",.libPaths()))
}

print('start')
x <- 2
while(x>3){
    print('A')
    seed <- 707
baseversion <- 'testmcmc1_'
nclusters <- 16L
niter <- 1024L*4L
niter0 <- 2048L
thin <- 4L
nstages <- 2L
ncheckprobs1 <- 16L
ncheckprobs2 <- 8L
maincov <- 'class'
family <- 'Palatino'
ndata <- 4096L
posterior <- TRUE
    x <- x-1
    print(x)
}
