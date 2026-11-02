CC = gcc
CFLAGS ?= -O2 -std=gnu11 -Wall -Wextra -fopenmp
LDLIBS = -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lcfitsio -lm

.PHONY: all clean
all: 1_ngrf/ngrf_2d 2b_visfitsgrid/visfitsgrid 3_vcg/vcg 4_fbe/viber 4_fbe/viber_equi_fast

1_ngrf/ngrf_2d: 1_ngrf/ngrf_2d.c 1_ngrf/ngrf_2d_funcs.c 1_ngrf/fitsprog_2d.c 1_ngrf/beam.c
	$(CC) $(CPPFLAGS) $(CFLAGS) 1_ngrf/ngrf_2d.c 1_ngrf/ngrf_2d_funcs.c 1_ngrf/fitsprog_2d.c 1_ngrf/beam.c $(LDFLAGS) $(LDLIBS) -o $@

2b_visfitsgrid/visfitsgrid: 2b_visfitsgrid/fitsprog.c 2b_visfitsgrid/visfitsgrid.c 2b_visfitsgrid/read_fits_func.c 2b_visfitsgrid/fitsprog.h 2b_visfitsgrid/read_fits_func.h
	$(CC) $(CPPFLAGS) $(CFLAGS) 2b_visfitsgrid/fitsprog.c 2b_visfitsgrid/visfitsgrid.c 2b_visfitsgrid/read_fits_func.c $(LDFLAGS) $(LDLIBS) -o $@

3_vcg/vcg: 3_vcg/read_fits_funcmfs.c 3_vcg/tgefuncs.c 3_vcg/main.c 3_vcg/tgefuncs.h 3_vcg/read_fits_func.h
	$(CC) $(CPPFLAGS) $(CFLAGS) 3_vcg/read_fits_funcmfs.c 3_vcg/tgefuncs.c 3_vcg/main.c $(LDFLAGS) $(LDLIBS) -o $@

4_fbe/viber: 4_fbe/viber.c 4_fbe/viber_funcs.c
	$(CC) $(CPPFLAGS) $(CFLAGS) 4_fbe/viber.c 4_fbe/viber_funcs.c $(LDFLAGS) $(LDLIBS) -o $@

4_fbe/viber_equi_fast: 4_fbe/viber_equi_fast.c 4_fbe/viber_funcs.c
	$(CC) $(CPPFLAGS) $(CFLAGS) 4_fbe/viber_equi_fast.c 4_fbe/viber_funcs.c $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f 1_ngrf/ngrf_2d 2b_visfitsgrid/visfitsgrid 3_vcg/vcg 4_fbe/viber 4_fbe/viber_equi_fast
