############################################################################
#                     							   						   #
#                      Basic operation mode of code                        #
#                                                                          #
############################################################################



OPT   +=   -DHBT  # Use Han2018 et al. HBT+ halo catalogue, close this default SUBFIND halo catalogue

#OPT   +=   -DUNBINDING  # Remove unbound particle use calculate spherical potential

#OPT   +=   -DOUTPE  # Use particle-particle method compute halo potential energy, when Np > 10^6 cost more time!!!!

#OPT   +=   -DHALOID  # Get halo particle id save in  halo_id_XXX.X  XXX is snap, Note: Unbinding can affect it, default long long type

#OPT   +=   -DPROJECTION  # Get around halo CIC projection(xy, yz, xz) save in halo_pic_XXX.X, cost more time!!!, close the UNBINDING

OPT   +=   -DDEBUG  # Debug every core output in Logfile ( Task_XXX, XXX is core number )

OPTIONS =   -DCOMPILETIMESETTINGS=\""$(OPT)"\"


SYSTYPE="ccc"


ifeq ($(SYSTYPE),"ccc")
CC       =   mpicc       
OPTIMIZE =   -O3 -Wall -g -m64  
GSL_INCL =  -I/home/tczhang/gsl-install/include
GSL_LIBS =  -L/home/tczhang/gsl-install/lib
FFTW_INCL=  -I/home/tczhang/fftw-install/include
FFTW_LIBS=  -L/home/tczhang/fftw-install/lib
MPICHLIB =
endif


EXEC   = CHOPPERS

OBJS   = main.o halo.o allvars.o peano.o density.o fitnfw.o unbinding.o spin.o shape.o vc.o energy.o projection.o

INCL   = allvars.h  proto.h  Makefile

LIBS   = $(GSL_LIBS) -lgsl -lgslcblas  -lm 

CFLAGS =    $(OPTIONS)  $(OPT) $(GSL_INCL) 

$(EXEC): $(OBJS) 
	$(CC) $(OPT) $(OBJS) $(LIBS)  -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)
