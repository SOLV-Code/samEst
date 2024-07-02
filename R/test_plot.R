data(harck)
library(samEst)

ricker_stan(harck)

f=samEst::ricker_stan(data=harck,type='static')

samEst::post_check(f,data=harck)
