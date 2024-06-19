#============================================
#Testing Thorson's 2024 AIc script
#implementing cAIC from Thorson 2024 paper
#https://pubmed.ncbi.nlm.nih.gov/38859712/
#Catarina wor
#August 2022
#============================================


devtools::document()
devtools::load_all()

#source("devs/calculate_EDF.R")
library(TMB)
#fit my TMB version
data(harck)
ptva<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1,priors_flag=0)
ptva$tmb_obj$fn()[1]
ptva$tmb_obj$report()$nll
ptva$tmb_obj$report()$renll
ptva$AICc
ptva$siga
ptva$sig


#thorson version
thdat<-list()
thdat$obs_S = harck$S
thdat$obs_logRS = harck$logRS
thdat$options_z<-c(1,0)

thparams <- list( "alpha" = 0,
                 "beta" = 0,
                 "ln_sigA" = 0,
                 "ln_sigB" = 0,
                 "ln_sigma" = 0,
                 "epsA_t" = rep(0,nrow(harck)),
                 "epsB_t" = rep(0,nrow(harck)) )
Map <- list()
Map$ln_sigB = factor(NA)
Map$epsB_t = factor( rep(NA,length(thparams$epsB_t)) )

TMB::compile('devs/Ricker_thorson.cpp')
dyn.load(dynlib('devs/Ricker_thorson'))

objth <- TMB::MakeADFun(data = thdat, 
      parameters = thparams, 
      map = Map,
      random = c( "epsA_t"), 
      DLL = "Ricker_thorson")

 
objth$env$beSilent()
nonvariance_fixed_effects = c("alpha","beta")

  # Optimize first time
optth <- nlminb( start=objth$par, obj=objth$fn, grad=objth$gr )
SD = sdreport( objth)
fit10 = list( "Opt"=optth,
              "SD" = SD,
              "Obj"=objth,
              "nonvariance_fixed_effects"=nonvariance_fixed_effects )

EDF10 = calculate_EDF( obj=fit10$Obj,
                      opt=fit10$Opt,
                      nonvariance_fixed_effects=fit10$nonvariance_fixed_effects,
                      refit = "full")

my_nonvariance_fixed_effects<-c("alphao","logbeta")
names(ptva)
myEDF10 = calculate_EDF( obj=ptva$tmb_obj,
                      opt=ptva$model,
                      prediction_name = "pred_logRS",
                      data_name = "obs_logRS",
                      nonvariance_fixed_effects=my_nonvariance_fixed_effects,
                      refit = "full")

p <- ricker_TMB(data=harck,priors_flag=0)
myEDF00 = calculate_EDF( obj=p$tmb_obj,
                      opt=p$model,
                      prediction_name = "pred_logRS",
                      data_name = "obs_logRS",
                      nonvariance_fixed_effects=c("alpha","logbeta"),
                      refit = "full")


length(p$model$par)

mymAIC = 2*ptva$model$obj + 2*length(ptva$model$par)
mycAIC = -2*sum(ptva$tmb_obj$report()$ll) + 2*myEDF10


cAIC10 = -2*sum(fit10$Obj$report()$loglik_t) + 2*EDF10

mAIC10 = 2*fit10$Opt$obj + 2*length(fit10$Opt$par)
mAICc10 = 2*fit10$Opt$obj + 2*length(fit10$Opt$par) +(2*length(fit10$Opt$par)*(length(fit10$Opt$par)+1)/(nrow(harck)-length(fit10$Opt$par)-1))

comp<-data.frame(my=c(mymAIC,mycAIC), th=c(mAIC10,cAIC10), type=c("marginal","cond"))
comp$diff=comp$my-comp$th


#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)
#library(samEst)

#read in data 
#use Harrison as an example


#sr <- read.csv("C:/Users/worc/Documents/timevarproject/simeval/data/samsimHarCk/HARSR.csv")

#head(sr)


#har<-data.frame(by=sr$Brood.Year,
#	S=sr$Sum.Total.Spawners,
#	R=sr$AEQ_Recruitment..age.2.5.,
#	logRS=log(sr$AEQ_Recruitment..age.2.5./sr$Sum.Total.Spawners))


#harck<-har[!is.na(har$S),]

#setwd("C:\\Users\\worc\\Documents\\timevarproject\\samEst\\data")
#usethis::use_data(harck)
#library(samEst)
data(harck)
plot(harck$S,harck$R)


##testing functions

p <- ricker_TMB(data=harck)
p$Smax
p$alpha

fullLL<-ricker_kf_TMB(data=harck,fullLL=TRUE)
reLL<-ricker_kf_TMB(data=harck)

fullLL$tmb_obj$fn()
reLL$tmb_obj$fn()

2*fullLL$tmb_obj$fn()+2*3

2*reLL$tmb_obj$fn()+2*3

pkf<-ricker_TMB(data=harck, AC=TRUE)


ip_logb_mean<-log(1/(max(harck$S)*.5))
ip_logb_sd<-sqrt(log(1+1))

p_ip <- ricker_TMB(data=harck,logb_p_mean=ip_logb_mean,logb_p_sd=ip_logb_sd)
p_ip$Smax
p_ip$alpha

simple_mod <- samEst::compile_code(type='static', ac=FALSE, par='n',lambertW = FALSE)
b <- ricker_stan(data=harck,iter = 800, mod=simple_mod)

pnp <- ricker_TMB(data=harck,prior=0)

pb <- ricker_stan(data=harck,iter = 2000)

mymod <- compile_code(type='static',ac=TRUE,par='n',caphigh=FALSE)
pb2 <- ricker_stan(data=harck,iter = 2000,AC=TRUE, mod = mymod)

names(pb)
p$alpha
pb$alpha

p$beta
pb$beta


pac<-ricker_TMB(data=harck, AC=TRUE)
pac[1:10]

lfostatic<-tmb_mod_lfo_cv(data=harck,model='static')
lfoac <- tmb_mod_lfo_cv(data=harck,model='staticAC')


sum(lfostatic)
sum(lfoac)
names(p)


ptva<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1)
ptva_ip<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1,logb_p_mean=ip_logb_mean,logb_p_sd=ip_logb_sd)

ptva[1:6]
ptva_ip[1:6]
pac$tmb_obj$report()$nll

pkfa<- ricker_kf_TMB(data=harck)
pkfa[1:6]

pkfa2<- ricker_kf_TMB(data=harck,fullLL=T)
pkfa2[1:6]

lfoalpha <- tmb_mod_lfo_cv(data=harck,tv.par='a', siglfo="obs")
sum(lfoalpha$lastparam)
sum(lfoalpha$last3paramavg)
sum(lfoalpha$last5paramavg)


phmm <- ricker_hmm_TMB(data=harck, tv.par='both')
phmm[1:10]

dirpr<-matrix(c(4,1,1,4),2,2)
phmm_dirpr <- ricker_hmm_TMB(data=harck, tv.par='both',dirichlet_prior=dirpr)
phmm_dirpr[1:10]


pbhmm <- ricker_hmm_stan(data, par='b')



phmmb <- ricker_hmm_TMB(data=harck, tv.par='b')
phmmb[1:8]


ptvb <- ricker_rw_TMB(data=harck,tv.par="b")

ptvab <- ricker_rw_TMB(data=harck,tv.par="both")



lfohmm <- tmb_mod_lfo_cv(data=harck,tv.par='HMM')

sum(lfohmm$regime_pick,na.rm=T)
sum(lfohmm$regime_average)


phmma <- ricker_hmm_TMB(data=harck, tv.par='a')

phmma <- ricker_hmm_TMB(data=harck, tv.par='a')


phmma[1:5]

lfohmma <- tmb_mod_lfo_cv(data=harck,tv.par='HMM_a')

sum(lfohmma$regime_pick)
sum(lfohmma$regime_average)

lfohmmb <- tmb_mod_lfo_cv(data=harck,tv.par='HMM_b')

sum(lfohmmb$regime_pick,na.rm=T)
sum(lfohmmb$regime_average)


phmmb <- ricker_hmm_TMB(data=harck, tv.par='b')
phmmb[1:5]




#stan functions

#lfo tmb testing
lfostatic<-tmb_mod_lfo_cv(data=harck,model='static', L=round((2/3)*nrow(harck)))
lfoac <- tmb_mod_lfo_cv(data=harck,model='staticAC', L=round((2/3)*nrow(harck)))
lfoalpha <- tmb_mod_lfo_cv(data=harck,model='rw_a', siglfo="obs", L=round((2/3)*nrow(harck)))
lfobeta <- tmb_mod_lfo_cv(data=harck,model='rw_b', siglfo="obs", L=round((2/3)*nrow(harck)))
lfoalphabeta <- tmb_mod_lfo_cv(data=harck,model='rw_both', siglfo="obs", L=round((2/3)*nrow(harck)))
lfohmma <- tmb_mod_lfo_cv(data=harck,model='HMM_a', L=round((2/3)*nrow(harck)))
lfohmmb <- tmb_mod_lfo_cv(data=harck,model='HMM_b', L=round((2/3)*nrow(harck)))
lfohmm <- tmb_mod_lfo_cv(data=harck,model='HMM', L=round((2/3)*nrow(harck)))


p <- ricker_TMB(data=harck)
pac<-ricker_TMB(data=harck, AC=TRUE)
ptva<- ricker_rw_TMB(data=harck,tv.par="a")
ptvb <- ricker_rw_TMB(data=harck,tv.par="b",sig_p_sd=1)
ptvab <- ricker_rw_TMB(data=harck,tv.par="both",sig_p_sd=.5)
phmma <- ricker_hmm_TMB(data=harck, tv.par='a')
phmmb <- ricker_hmm_TMB(data=harck, tv.par='b')
phmm <- ricker_hmm_TMB(data=harck, tv.par='both')




c(
p$AICc,
pac$AICc,
ptva$AICc,
ptvb$AICc, 
ptvab$AICc,
phmma$AICc,
phmmb$AICc,
phmm$AICc)

c(
p$BIC,
pac$BIC,
ptva$BIC,
ptvb$BIC, 
ptvab$BIC,
phmma$BIC,
phmmb$BIC,
phmm$BIC)




#=====================================================
#RTMB


devtools::document()
devtools::load_all()
#library(RTMB)


p <- ricker_TMB(data=harck)
p2 <- ricker_RTMB(data=harck)

p$alpha
p2$logalpha

pac <- ricker_TMB(data=harck,AC=TRUE,priors_flag=0)
pac$rho
pac$alpha
pac$beta
pac$Smax
pac$tmb_obj$fn()
names(pac)

devtools::document()
devtools::load_all()
pac2 <- ricker_RTMB(data=harck,AC=TRUE,priors_flag=0)
pac2$rho
pac2$logalpha
pac2$beta
pac2$Smax
pac2$sigAR

pac2$obj$fn()


data=harck
dat<- list(
    obs_S = as.numeric(data$S),
    obs_logRS = data$logRS,
    by = data$by,
    priors_flag=1,
    stan_flag=0,
    sig_p_sd=1,
    logb_p_mean=-12,
    logb_p_sd=3
  )
  

  magS <- log10_ceiling(max(data$S))
  initlm <- lm(obs_logRS~obs_S, data=dat)

 
param <- list(
    logalpha   = initlm$coefficients[[1]],
    logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
    logsigobs = log(1),
    ar1_phi=3.0)


library(RTMB)
ricker_ac_RTMB_fn(param,dat)


tmbfn<-function(param){ricker_ac_RTMB_fn(param,dat=dat)}
obj <- MakeADFun(tmbfn, param, silent=TRUE)


opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)


p$alpha
obj$report()$logalpha

p$Smax
obj$report()$Smax



 final_grads <- sdr$gradient.fixed
  bad_eig <- FALSE
  conv_problem<-FALSE
  if (!is.null(sdr$pdHess)) {
    if (!sdr$pdHess) {
      warning("The model may not have converged: ",
        "non-positive-definite Hessian matrix.", call. = FALSE)
      conv_problem<-TRUE
    } else {
      eigval <- try(1 / eigen(sdr$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
          "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
        conv_problem<-TRUE
      }
      if (any(final_grads > 0.01)){
        warning("The model may not have converged. ",
          "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
        conv_problem<-TRUE



(Smsy <- (1 - gsl::lambert_W0(exp(1 - p$alpha))) /p$beta)
(Smsy2 <- (1 - samest_lambertW(exp(1 - p$alpha))) /p$beta)




pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")


