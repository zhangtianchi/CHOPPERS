
CC       =  mpicc               # sets the C-compiler

OPTIMIZE =  -O2 -Wall -g   # sets optimization and warning flags

OPTIONS =  $(OPTIMIZE) $(OPT)

EXEC   = HBT2LSubfind

OBJS   = main.o 

INCL   =  Makefile


GSL_INCL =  -I/home/tczhang/gsl-install/include
GSL_LIBS =  -L/home/tczhang/gsl-install/lib

CFLAGS = $(OPTIONS) $(GSL_INCL)


LIBS   =  -lgsl -lgslcblas -lm  $(GSL_LIBS)

# Implicit Rules
.SUFFIXES: .o .c


.c.o:
	$(CC) -c $< -o $*.o $(CFLAGS)


$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f *~ $(OBJS) $(EXEC)

