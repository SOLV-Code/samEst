all: Ricker_simple.so 

Ricker_simple.so: Ricker_simple.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_simple.cpp')"

clean:
	rm -rf *o
