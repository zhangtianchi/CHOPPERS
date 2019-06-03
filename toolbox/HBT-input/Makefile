





EXEC   = HBT2HaloTrace

GSL_INCL =  -I/home/tczhang/gsl-install/include
GSL_LIBS =  -L/home/tczhang/gsl-install/lib

HDF5INCL =  -I/data/dell1/userdir/tczhang/hdf5/include
HDF5LIB  =  -L/data/dell1/userdir/tczhang/hdf5/lib  -lhdf5



OBJS   = main.o readhdf.o allvars.o 


INCL   = allvars.h   proto.h Makefile



OPTIONS =  -O3 -Wall -m64  #${GSL_INCL}	${HDF5INCL}
CFLAGS = $(OPTIONS) $(GSL_INCL)  $(HDF5INCL)




LIBS = $(HDF5LIB) -g    $(GSL_LIBS) -lgsl -lgslcblas -lm 


#CC     =  gcc


$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)
