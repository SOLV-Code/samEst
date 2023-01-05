all: Ricker_simple.so Ricker_autocorr.so Ricker_tva.so Ricker_tvlogb.so Ricker_tva_tvb.so SR_HMM.so SR_HMM_a.so SR_HMM_b.so Rickerkf_fullLL.so Rickerkf.so

Ricker_simple.so: Ricker_simple.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_simple.cpp')"

Ricker_autocorr.so: Ricker_autocorr.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_autocorr.cpp')"

Ricker_tva.so: Ricker_tva.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_tva.cpp')"

Ricker_tvlogb.so: Ricker_tvlogb.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_tvlogb.cpp')"

Ricker_tva_tvb.so: Ricker_tva_tvb.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_tva_tvb.cpp')"

SR_HMM.so: SR_HMM.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM.cpp')"

SR_HMM_a.so: SR_HMM_a.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_a.cpp')"

SR_HMM_b.so: SR_HMM_b.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_b.cpp')"


Rickerkf_fullLL.so: Rickerkf_fullLL.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf_fullLL.cpp')"


Rickerkf.so: Rickerkf.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf.cpp')"



clean:
	rm -rf *o
