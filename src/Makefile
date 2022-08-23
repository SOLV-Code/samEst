all: Ricker_simple.so Ricker_tva.so 

Ricker_simple.so: Ricker_simple.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_simple.cpp')"

Ricker_tva.so: Ricker_tva.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_tva.cpp')"

clean:
	rm -rf *o
