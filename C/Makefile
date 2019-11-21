ifdef CFITSIO_PATH
INCLUDES += -I $(CFITSIO_PATH)/include
LIBS += -L $(CFITSIO_PATH)/lib -Xlinker -rpath -Xlinker $(CFITSIO_PATH)
endif
LIBS += -lcfitsio

ifdef NETCDF_PATH
INCLUDES += -I $(NETCDF_PATH)/include
LIBS += -L $(NETCDF_PATH)/lib -Xlinker -rpath -Xlinker $(NETCDF_PATH)
endif
LIBS += -lnetcdf

LIBS += -lm

all:	spec_driver_fits

clean:
	rm -f *.o *~ spec_driver_fits

Cube.o : Cube.h Cube.c
	gcc -c Cube.c -O3 $(INCLUDES)

Plane.o : Plane.h Plane.c
	gcc -c Plane.c -O3 $(INCLUDES)

SpecFile.o : SpecFile.h SpecFile.c
	gcc -c SpecFile.c -O3 $(INCLUDES)

ConvolveFunction.o : ConvolveFunction.h ConvolveFunction.c
	gcc -c ConvolveFunction.c -O3

OTFParameters.o : OTFParameters.h OTFParameters.c
	gcc -c OTFParameters.c -O3

spec_driver_fits : spec_driver_fits.c Cube.o Plane.o SpecFile.o ConvolveFunction.o OTFParameters.o
	gcc -o spec_driver_fits -O3 spec_driver_fits.c Cube.o Plane.o SpecFile.o ConvolveFunction.o OTFParameters.o $(INCLUDES) $(LIBS) 

