all: EndTB.so

EndTB.so: EndTB.cpp
	Rscript --vanilla -e "TMB::compile('EndTB.cpp', '-I../inst/include')"

clean:
	rm -rf *.o