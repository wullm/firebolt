#Compiler options
GCC = gcc

#Libraries
INI_PARSER = parser/minIni.o
STD_LIBRARIES = -lm
FFTW_LIBRARIES = -lfftw3
HDF5_LIBRARIES = -lhdf5
GSL_LIBRARIES = -lgsl -lgslcblas

GSL_INCLUDES =

HDF5_INCLUDES += -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
HDF5_LIBRARIES += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -I/usr/include/hdf5/serial

#Putting it together
INCLUDES = $(HDF5_INCLUDES) $(GSL_INCLUDES)
LIBRARIES = $(INI_PARSER) $(STD_LIBRARIES) $(FFTW_LIBRARIES) $(HDF5_LIBRARIES) $(GSL_LIBRARIES)
CFLAGS = -Wall -Ofast -march=native -fopenmp -fPIC -ggdb3

OBJECTS = lib/*.o

all:
	make minIni
	$(GCC) src/input.c -c -o lib/input.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/output.c -c -o lib/output.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/background.c -c -o lib/background.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/background_interp.c -c -o lib/background_interp.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/perturb_data.c -c -o lib/perturb_data.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/perturb_interp.c -c -o lib/perturb_interp.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/ic.c -c -o lib/ic.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/evolve.c -c -o lib/evolve.o $(INCLUDES) $(CFLAGS)

	$(GCC) src/fft.c -c -o lib/fft.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/multipoles.c -c -o lib/multipoles.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/multipole_interp.c -c -o lib/multipole_interp.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/grids.c -c -o lib/grids.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/evaluate.c -c -o lib/evaluate.o $(INCLUDES) $(CFLAGS)

	$(GCC) src/firebolt.c -o firebolt $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)
	$(GCC) src/firebolt_single.c -o firebolt_single $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)
	$(GCC) src/firebolt_print.c -o firebolt_print $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)
	$(GCC) src/firebolt_check.c -o firebolt_check $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)
	$(GCC) src/firebolt_render.c -o firebolt_render $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)

	$(GCC) -shared -o libfirebolt.so $(INCLUDES) lib/multipoles.o lib/multipole_interp.o lib/evolve.o lib/fft.o lib/grids.o lib/evaluate.o lib/input.o lib/perturb_data.o lib/perturb_interp.o $(LIBRARIES) $(CFLAGS)

	$(GCC) src/testlib.c -o testlib -L/home/qvgd89/firebolt -lfirebolt -lfftw3 -lgsl -lgslcblas $(HDF5_LIBRARIES) $(INCLUDES) $(CFLAGS) -lm -Wl,-rpath=/home/qvgd89/firebolt

minIni:
	cd parser && make

check:
	cd tests && make
