#============================================
#commands used for packages development and testing
#Catarina wor
#August 2022
#============================================



devtools::document()
devtools::load_all()





#read in data 
#use Harrison as an example


sr<-read.csv("C:/Users/worc/Documents/timevarproject/simeval/data/samsimHarCk/HARSR.csv")

head(sr)


har<-data.frame(by=sr$Brood.Year,
	S=sr$Sum.Total.Spawners,
	R=sr$AEQ_Recruitment..age.2.5.,
	logRS=log(sr$AEQ_Recruitment..age.2.5./sr$Sum.Total.Spawners))


#harck<-har[!is.na(har$S),]

#setwd("C:\\Users\\worc\\Documents\\timevarproject\\samEst\\data")
#usethis::use_data(harck)

data(harck)
plot(harck$S,harck$R)


##testing functions

p<-rickerTMB(data=harck)
pb<-rickerstan(data=harck,iter = 2000)


p$alpha
pb$alpha

p$beta
pb$beta


pac<-rickerTMB(data=harck, AC=TRUE)

lfostatic<-tmb_mod_lfo_cv(data=harck,tv.par='static')
lfoac <- tmb_mod_lfo_cv(data=harck,tv.par='staticAC')


sum(lfostatic)
sum(lfoac)
names(p)

ptva <- ricker_rwa_TMB(data=harck)
ptva2 <- ricker_rw_TMB(data=harck,tvpar="a")

plot(ptva$alpha, type="b", lwd=2)
lines(ptva2$alpha, type="b",lty=2, lwd=2, col="red")


lfoalpha <- tmb_mod_lfo_cv(data=harck,tv.par='alpha', siglfo="obs")
sum(lfoalpha$lastparam)
sum(lfoalpha$last3paramavg)
sum(lfoalpha$last5paramavg)


phmm <- ricker_HMM_TMB(data=harck)

phmm[1:10]
fit_past_hmm_tmb[1:7]

lfohmm <- tmb_mod_lfo_cv(data=harck,tv.par='HMM')

sum(lfohmm$regime_pick,na.rm=T)
sum(lfohmm$regime_average)


ptvb <- ricker_rwb_TMB(data=harck)
ptvb2 <- ricker_rw_TMB(data=harck,tvpar="b")


plot(ptvb$beta, type="b", lwd=2)
lines(ptvb2$beta, type="b",lty=2, lwd=2, col="red")


ptvab <- ricker_rwab_TMB(data=harck)
ptvab2 <- ricker_rw_TMB(data=harck,tvpar="both")


plot(ptvab$beta, type="b", lwd=2)
lines(ptvab2$beta, type="b",lty=2, lwd=2, col="red")

plot(ptvab$alpha, type="b", lwd=2)
lines(ptvab2$alpha, type="b",lty=2, lwd=2, col="red")


phmma<-ricker_HMM_TMB_a(data=harck)
phmma[1:7]
lfohmma <- tmb_mod_lfo_cv(data=harck,tv.par='HMM_a')

sum(lfohmma$regime_pick)
sum(lfohmma$regime_average)

lfohmmb <- tmb_mod_lfo_cv(data=harck,tv.par='HMM_b')

sum(lfohmmb$regime_pick,na.rm=T)
sum(lfohmmb$regime_average)



phmmb<-ricker_HMM_TMB_b(data=harck)
phmmb[1:7]
#stan functions








