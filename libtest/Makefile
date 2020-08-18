#Compiler options
GCC = gcc

#Libraries
INI_PARSER = parser/minIni.o
STD_LIBRARIES = -lm
FFTW_LIBRARIES = -lfftw3
HDF5_LIBRARIES = -lhdf5
GSL_LIBRARIES = -lgsl -lgslcblas

GSL_INCLUDES =

FIREBOLT_PATH = "../"

HDF5_INCLUDES += -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
HDF5_LIBRARIES += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -I/usr/include/hdf5/serial

#Putting it together
INCLUDES = $(HDF5_INCLUDES) $(GSL_INCLUDES)
LIBRARIES = $(INI_PARSER) $(STD_LIBRARIES) $(FFTW_LIBRARIES) $(HDF5_LIBRARIES) $(GSL_LIBRARIES)
CFLAGS = -Wall -Ofast -march=native -fopenmp -fPIC -ggdb3

OBJECTS = ../lib/*.o obj/*.o

all:
	make minIni

	$(GCC) src/input.c -c -o obj/input.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/output.c -c -o obj/output.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/background.c -c -o obj/background.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/background_interp.c -c -o obj/background_interp.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/perturb_data.c -c -o obj/perturb_data.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/perturb_interp.c -c -o obj/perturb_interp.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/ic.c -c -o obj/ic.o $(INCLUDES) $(CFLAGS)

	$(GCC) src/firebolt_check.c -o firebolt_check -I/home/qvgd89/firebolt/include $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS)

	$(GCC) src/testlib.c -o testlib obj/perturb_data.o obj/perturb_interp.o obj/input.o -L$(FIREBOLT_PATH) -I$(FIREBOLT_PATH)/include -lfirebolt -lfftw3 -lgsl -lgslcblas $(INI_PARSER) $(HDF5_LIBRARIES) $(INCLUDES) $(CFLAGS) -lm -Wl,-rpath=$(FIREBOLT_PATH)

minIni:
	cd parser && make

check:
	cd tests && make
