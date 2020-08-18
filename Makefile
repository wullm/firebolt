#Compiler options
GCC = gcc

#Libraries
STD_LIBRARIES = -lm
FFTW_LIBRARIES = -lfftw3
HDF5_LIBRARIES = -lhdf5
GSL_LIBRARIES = -lgsl -lgslcblas

GSL_INCLUDES =

HDF5_INCLUDES += -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
HDF5_LIBRARIES += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -I/usr/include/hdf5/serial

#Putting it together
INCLUDES = $(HDF5_INCLUDES) $(GSL_INCLUDES)
LIBRARIES = $(STD_LIBRARIES) $(FFTW_LIBRARIES) $(HDF5_LIBRARIES) $(GSL_LIBRARIES)
CFLAGS = -Wall -Wshadow -Ofast -march=native -fopenmp -fPIC -ggdb3

OBJECTS = lib/*.o

all:
	$(GCC) src/evolve.c -c -o lib/evolve.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/fft.c -c -o lib/fft.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/multipoles.c -c -o lib/multipoles.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/multipole_interp.c -c -o lib/multipole_interp.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/grids.c -c -o lib/grids.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/evaluate.c -c -o lib/evaluate.o $(INCLUDES) $(CFLAGS)

	$(GCC) -shared -o libfirebolt.so $(INCLUDES) lib/multipoles.o lib/multipole_interp.o lib/evolve.o lib/fft.o lib/grids.o lib/evaluate.o $(LIBRARIES) $(CFLAGS)

	# $(GCC) src/testlib.c -o testlib lib/perturb_data.o lib/perturb_interp.o lib/input.o -L/home/qvgd89/firebolt -lfirebolt -lfftw3 -lgsl -lgslcblas $(HDF5_LIBRARIES) $(INCLUDES) $(CFLAGS) -lm -Wl,-rpath=/home/qvgd89/firebolt

check:
	cd libtest && make
