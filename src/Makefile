all: Ricker_simple.so Ricker_autocorr.so SR_HMM.so SR_HMM_a.so SR_HMM_b.so Rickerkf.so Ricker_tv_all.so

Ricker_simple.so: Ricker_simple.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_simple.cpp')"

Ricker_autocorr.so: Ricker_autocorr.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_autocorr.cpp')"
		
SR_HMM.so: SR_HMM.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM.cpp')"

SR_HMM_a.so: SR_HMM_a.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_a.cpp')"

SR_HMM_b.so: SR_HMM_b.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_b.cpp')"

Rickerkf.so: Rickerkf.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf.cpp')"

Ricker_tv_all.so: Ricker_tv_all.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_tv_all.cpp')"


clean:
	rm -rf *o
