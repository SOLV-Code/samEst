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


p<-rickerTMB(data=harck)

ptva <- ricker_rwa_TMB(data=harck)

plot(ptva$alpha, type="b")

phmm <- ricker_HMM_TMB(data=harck)


ptvb <- ricker_rwb_TMB(data=harck)