all: Ricker_simple.so Ricker_tva.so Ricker_tvlogb.so Ricker_tva_tvb.so SR_HMM.cpp SR_HMM_a.cpp SR_HMM_b.cpp

Ricker_simple.so: Ricker_simple.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_simple.cpp')"

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

clean:
	rm -rf *o
