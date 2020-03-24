#Compiler options
GCC = gcc

#Libraries
INI_PARSER = ../parser/minIni.o
STD_LIBRARIES = -lm
HDF5_LIBRARIES = -lhdf5
GSL_LIBRARIES = -lgsl -lgslcblas

GSL_INCLUDES =

HDF5_INCLUDES += -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
HDF5_LIBRARIES += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -I/usr/include/hdf5/serial

#Putting it together
INCLUDES = $(HDF5_INCLUDES) $(GSL_INCLUDES)
LIBRARIES = $(INI_PARSER) $(STD_LIBRARIES) $(HDF5_LIBRARIES) $(GSL_LIBRARIES)
CFLAGS = -Wall

OBJECTS = ../lib/*.o

all:
	@#$(GCC) test_minIni.c -o test_minIni $(INI_PARSER)
	@#@./test_minIni

	$(GCC) test_input.c -o test_input $(OBJECTS) $(LIBRARIES) $(CFLAGS) $(INCLUDES)
	@./test_input

	$(GCC) test_background.c -o test_background $(OBJECTS) $(LIBRARIES) $(CFLAGS) $(INCLUDES)
	@./test_background

	$(GCC) test_perturb_data.c -o test_perturb_data $(OBJECTS) $(LIBRARIES) $(CFLAGS) $(INCLUDES)
	@./test_perturb_data
