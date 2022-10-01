#============================================
#commands used for packages development and testing
#Catarina wor
#August 2022
#============================================



devtools::document()
devtools::load_all()





#read in data 
#use Harrison as an example


sr <- read.csv("C:/Users/worc/Documents/timevarproject/simeval/data/samsimHarCk/HARSR.csv")

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

p <- ricker_TMB(data=harck)

pb <- ricker_stan(data=harck,iter = 2000)


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


ptva<- ricker_rw_TMB(data=harck,tv.par="a")
ptva[1:10]

lfoalpha <- tmb_mod_lfo_cv(data=harck,tv.par='a', siglfo="obs")
sum(lfoalpha$lastparam)
sum(lfoalpha$last3paramavg)
sum(lfoalpha$last5paramavg)


phmm <- ricker_hmm_TMB(data=harck, tv.par='both')
phmm[1:8]

pbhmm <- ricker_hmm_stan(data, par='both')

pb <- rickerstan(data=harck,iter = 2000)







ptvb <- ricker_rw_TMB(data=harck,tv.par="b")



ptvab <- ricker_rw_TMB(data=harck,tv.par="both")



lfohmm <- tmb_mod_lfo_cv(data=harck,tv.par='HMM')

sum(lfohmm$regime_pick,na.rm=T)
sum(lfohmm$regime_average)


phmma <- ricker_HMM_TMB(data=harck, tv.par='a')
phmma[1:5]

phmma2<-ricker_HMMa_TMB(data=harck)
phmma2[1:5]
lfohmma <- tmb_mod_lfo_cv(data=harck,tv.par='HMM_a')

sum(lfohmma$regime_pick)
sum(lfohmma$regime_average)

lfohmmb <- tmb_mod_lfo_cv(data=harck,tv.par='HMM_b')

sum(lfohmmb$regime_pick,na.rm=T)
sum(lfohmmb$regime_average)


phmmb <- ricker_HMM_TMB(data=harck, tv.par='b')
phmmb[1:5]


phmmb2<-ricker_HMMb_TMB(data=harck)
phmmb2[1:5]
#stan functions








