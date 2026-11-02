FITSDIR=/usr/local
target: visfitsgrid clean tidy visfitsgrid_sampling visfitsgrid_BAIC

visfitsgrid: fitsprog.c visfitsgrid.c read_fits_func.c 
	gcc -g -o visfitsgrid fitsprog.c visfitsgrid.c read_fits_func.c -I$(FITSDIR)/include -L$(FITSDIR)/lib -lcfitsio -lfftw3 -lm
	rm -rf *~
	

clean:
	rm -rf visfitsgrid visfitsgrid_sampling *.FITS

tidy:
	rm -rf *~
	rm -rf *.o

