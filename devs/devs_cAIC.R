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
p <- ricker_TMB(data=harck,sig_p_sd=1,priors_flag=0)
ptva<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1,priors_flag=0)
ptvb<- ricker_rw_TMB(data=harck,tv.par="b",sig_p_sd=1,priors_flag=0)

ptvap<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1,priors_flag=1)
ptvbp<- ricker_rw_TMB(data=harck,tv.par="b",sig_p_sd=1,priors_flag=1)

ptvam<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1,priors_flag=1,
 AICc_type="marginal")

p$AICc
ptva$AICc
ptva$EDF

ptvb$AICc
ptvap$AICc
ptvbp$AICc
ptvam$AICc



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


#======================================


#test with simulated data

simdat<-readRDS("C:/Users/worc/Documents/timevar/simest-tv/outs/SamSimOutputs/simData/stationary/stationary/stationary_fixedER_CUsrDat.RData")$srDatout

u=1

dat <- simdat[simdat$iteration==u,]
dat <- dat[dat$year>(max(dat$year)-46),]
dat <- dat[!is.na(dat$obsRecruits),]
df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

Smax_mean<-(max(df$S)*.5)
Smax_sd<-Smax_mean
 
logbeta_pr_sig=sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
logbeta_pr=log(1/(Smax_mean))-0.5*logbeta_pr_sig^2


#options
silent = FALSE 
control = TMBcontrol()
ini_param=NULL
tmb_map = list()
priors_flag=1
stan_flag=0
sig_p_sd=1
siga_p_sd=1
sigb_p_sd=.3
logb_p_mean=logbeta_pr
logb_p_sd=logbeta_pr_sig
  
AICc_type=c("conditional", "marginal")[1]
deltaEDF=0.0001

tmb_data <- list(
    obs_S = df$S,
    obs_logRS = df$logRS,
    priors_flag=priors_flag,
    stan_flag=stan_flag,
    sig_p_sd=sig_p_sd,
    logb_p_mean=logb_p_mean,
    logb_p_sd=logb_p_sd
  )

  
magS <- log10_ceiling(max(df$S))
initlm<-lm(logRS~S, data=df)
  
tmb_data$siga_p_sd=siga_p_sd


tmb_params <- list(alphao   = initlm$coefficients[[1]],
                   logbeta = ifelse(initlm$coefficients[[2]]>0,
                                   log(1/magS),
                                   log(-initlm$coefficients[[2]])),
                   logsigobs = log(.5),
                   logsiga = log(.5),
                   alpha = rep(initlm$coefficients[[1]],length(tmb_data$obs_S)))
    
tmb_random <- "alpha"
tmb_obj <- TMB::MakeADFun(data = tmb_data, 
                              parameters = tmb_params, 
                              map = tmb_map,
                              random = tmb_random, 
                              DLL = "Ricker_tva", 
                              silent = silent)

lowlimit <- c(0.01,-20,log(0.01),log(0.01))
hightlimit <- c(20,-4,log(2),log(2))
    
nonvariance_fixed_effects<-c("alphao","logbeta")

clss <- "Ricker_tva"
npar <- 4
npar_all <- 4+(length(df$S)-1)

  
 

#===================================
# TMB fit
#===================================

tmb_opt <- stats::nlminb(
    start = tmb_obj$par, 
    objective = tmb_obj$fn, 
    gradient = tmb_obj$gr,
    control = control,
    lower = lowlimit, 
    upper = hightlimit)

grads_nstp<-matrix(NA, ncol=4,nrow=20)
myEDF<- rep(NA,20)
mypars<-matrix(NA, ncol=4,nrow=20)
  
sd_report <- TMB::sdreport(tmb_obj)
conv <- get_convergence_diagnostics(sd_report)
grads_nstp[1,]<-conv$final_grads
mypars[1,] <-tmb_opt$par
fittoedf = list( "Opt"=tmb_opt,
              "SD" = sd_report,
              "Obj"=tmb_obj,
              "nonvariance_fixed_effects"=nonvariance_fixed_effects )

  
myEDF[1] = calculate_EDF( obj=fittoedf$Obj,
                      opt=fittoedf$Opt,
                      prediction_name = "pred_logRS",
                      data_name = "obs_logRS",
                      delta = deltaEDF,
                      nonvariance_fixed_effects=fittoedf$nonvariance_fixed_effects,
                      refit = "full")



#newton steps
for(i in seq_len(19)){
  g = as.numeric( tmb_obj$gr(tmb_opt$par) )
  h = optimHess(tmb_opt$par, fn=tmb_obj$fn, gr=tmb_obj$gr)
  tmb_opt$par = tmb_opt$par - solve(h, g)
  tmb_opt$objective = tmb_obj$fn(tmb_opt$par)
  sd_report <- TMB::sdreport(tmb_obj)
  conv <- get_convergence_diagnostics(sd_report)
  grads_nstp[i+1,]<-conv$final_grads
  mypars[i+1,] <-tmb_opt$par
  fittoedf = list( "Opt"=tmb_opt,
              "SD" = sd_report,
              "Obj"=tmb_obj,
              "nonvariance_fixed_effects"=nonvariance_fixed_effects )
  
  
  myEDF[i+1] = calculate_EDF( obj=fittoedf$Obj,
                      opt=fittoedf$Opt,
                      prediction_name = "pred_logRS",
                      data_name = "obs_logRS",
                      delta = deltaEDF,
                      nonvariance_fixed_effects=fittoedf$nonvariance_fixed_effects,
                      refit = "full")
}



fittoedf = list( "Opt"=tmb_opt,
              "SD" = sd_report,
              "Obj"=tmb_obj,
              "nonvariance_fixed_effects"=nonvariance_fixed_effects )

  
myEDF = calculate_EDF( obj=fittoedf$Obj,
                      opt=fittoedf$Opt,
                      prediction_name = "pred_logRS",
                      data_name = "obs_logRS",
                      delta = deltaEDF,
                      nonvariance_fixed_effects=fittoedf$nonvariance_fixed_effects,
                      refit = "full")



myEDF2 = calculate_EDF( obj=tmb_obj,
                      opt=tmb_opt,
                      prediction_name = "pred_logRS",
                      data_name = "obs_logRS",
                      delta = deltaEDF,
                      nonvariance_fixed_effects=fittoedf$nonvariance_fixed_effects,
                      refit = "full")
    













palla0<-ricker_rw_TMBall(data=df,tv.par='a',priors_flag=0) 
palla0$alpha
palla0$tmb_obj$fn()[1]

ptva0 <- ricker_rw_TMB(data=df,tv.par='a',priors_flag=0)
ptva0$alpha
ptva0$tmb_obj$fn()[1]


 
p0 <- ricker_TMB(data=df,priors_flag=0)
ptvb0 <- ricker_rw_TMB(data=df,tv.par='b',priors_flag=0)


p0$AICc
ptva0$AICc
ptvb0$AICc
 
ptva0$EDF

p <- ricker_TMB(data=df,logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
ptva <- ricker_rw_TMB(data=df,tv.par='a',logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
ptvb <- ricker_rw_TMB(data=df, tv.par='b',sigb_p_sd=1,logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)

ptvb$model$convergence
ptvb$conv_problem
 


ptvb0$model$convergence
ptvb0$conv_problem
 

p$AICc
ptva$AICc
ptvb$AICc









