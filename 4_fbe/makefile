target:	viber viber_equi_fast PS_2D_VCG tidy clean

CC=gcc
CFLAGS=-g

########## change the following paths according to your setup ##########
FFTW_LIB= -L/usr/local/lib/
CFIO_LIB= -L/usr/local/lib/
GSL_LIB= -L/usr/local/lib/

FFTW_INCLUDE= -I/usr/local/include/
CFIO_INCLUDE= -I/usr/local/include/
GSL_INCLUDE= -I/usr/local/include/
#####################################################


viber:	viber.c viber_funcs.c
	$(CC) $(CFLAGS) $(FFTW_INCLUDE) $(CFIO_INCLUDE) $(GSL_INCLUDE) viber.c viber_funcs.c -fopenmp -lm $(FFTW_LIB) -lfftw3 -lfftw3_omp $(CFIO_LIB) -lcfitsio $(GSL_LIB) -lgsl -lgslcblas -o viber
	rm -rf *~ *.o
	

viber_equi_fast:	viber_equi_fast.c viber_funcs.c
	$(CC) $(CFLAGS) $(FFTW_INCLUDE) $(CFIO_INCLUDE) $(GSL_INCLUDE) viber_equi_fast.c viber_funcs.c -fopenmp -lm $(FFTW_LIB) -lfftw3 -lfftw3_omp $(CFIO_LIB) -lcfitsio $(GSL_LIB) -lgsl -lgslcblas -o viber_equi_fast
	rm -rf *~ *.o
	
        
clean:
	rm -rf viber viber_equi_fast PS_2D_VCG masked_* *~ *.o 
        
tidy:
	rm -rf masked_* modes *~ *.o bsout* apsout*
	
	
