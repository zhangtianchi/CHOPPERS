############################################################################
#                     							   						   #
#                      Basic operation mode of code                        #
#                                                                          #
############################################################################



OPT   +=   -DHBT         # Use Han et al 2018. HBT+ halo catalogue, close this default SUBFIND halo catalogue

#OPT   +=   -DUNBINDING  # Remove unbound particle use calculate spherical potential

#OPT   +=   -DOCTREE      # Use octree method calculate potential, otherwise use particle-particle method compute halo potential, random pick 1000 particles calculate potential, then weighted ( see Neto et. al 2007  )

#OPT   +=   -DHALOID     # Get halo particle id save in  halo_id_XXX.X  XXX is snap, Note: Unbinding can affect it, default long long type

#OPT   +=   -DPROJECTION  # Get around halo CIC projection(xy, yz, xz) save in halo_pic_XXX.X, cost more time!!!, please close the UNBINDING first

#OPT   +=   -DDEBUG  # Debug every core output in Logfile ( Task_XXX, XXX is core number )

#OPT   +=   -DSIMPLIFY       # simplify version only calculate mass (drop velocity properties, concenreation, density), very fast!

#OPT   +=   -DDeltaVir       #  use  mvir rvir,  by  Bryan1998 fit equation replace r200 and m200.

#OPT   +=   -DUSEHBTR200     #  when halo become subhalo, if use 200rhoc, mass become very large, we use HBT2 orgin r200

#OPT   +=   -DKERNELDP       #  if you want Calc KernelDP( Reed2005 AppendixA Method  ) Open it!
########################################################################################################################

#OPT   +=   -DMYWORK      #  This is for my work (output density profile, range: [softening,r200], 20 logbins, physical unit, Note: Must turn off unbinding, only concenreation is right, every properties is not accuracy) , please turn off the option if you use this code.





OPTIONS =   -DCOMPILETIMESETTINGS=\""$(OPT)"\"


SYSTYPE="nova"


ifeq ($(SYSTYPE),"nova")
CC       =   mpicc       
OPTIMIZE =   -O3 -Wall -g -m64  
GSL_INCL =  -I/home/tczhang/gsl-install/include
GSL_LIBS =  -L/home/tczhang/gsl-install/lib
FFTW_INCL=  -I/home/tczhang/fftw-install/include
FFTW_LIBS=  -L/home/tczhang/fftw-install/lib
MPICHLIB =
endif

#CC       =   mpicc        # sets the C-compiler (default)
#OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)

EXEC   =  CHOPPERS

OBJS   = main.o halo.o allvars.o peano.o density.o fitnfw.o unbinding.o spin.o shape.o vc.o energy.o projection.o tree.o

INCL   = allvars.h  proto.h  Makefile

LIBS   = $(GSL_LIBS) -lgsl -lgslcblas  -lm 

CFLAGS =    $(OPTIONS)  $(OPT) $(GSL_INCL) 

$(EXEC): $(OBJS) 
	$(CC) $(OPT) $(OBJS) $(LIBS)  -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)
