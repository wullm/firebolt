#Compiler options
GCC = gcc

#Libraries
INI_PARSER = parser/minIni.o
STD_LIBRARIES = -lm
HDF5_LIBRARIES = -lhdf5
GSL_LIBRARIES = -lgsl -lgslcblas

GSL_INCLUDES =

HDF5_INCLUDES += -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
HDF5_LIBRARIES += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -I/usr/include/hdf5/serial

#Putting it together
INCLUDES = $(HDF5_INCLUDES) $(GSL_INCLUDES)
LIBRARIES = $(INI_PARSER) $(STD_LIBRARIES) $(HDF5_LIBRARIES) $(GSL_LIBRARIES)
CFLAGS = -Wall -Ofast -march=native -fopenmp

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
	$(GCC) src/firebolt.c -o firebolt $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)
	$(GCC) src/firebolt_single.c -o firebolt_single $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)
	$(GCC) src/firebolt_print.c -o firebolt_print $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)
	$(GCC) src/firebolt_check.c -o firebolt_check $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)

minIni:
	cd parser && make

check:
	cd tests && make
