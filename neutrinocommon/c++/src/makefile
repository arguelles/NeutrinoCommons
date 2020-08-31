# Compiler
GCC=g++-4.5

# Directories
LIBDIR        = ../lib
INCDIR        = ../inc

CFLAGS        = -Wall -g -O3 
INCFLAGS      = -I$(INCDIR)
LIBFLAGS      = 
LDFLAGS       = -lm -lgsl -lgslcblas -lstdc++ -lrt -lpthread

# Project files
LIBNAME = libneutrinocommon.a
OBJ     = body.o neuosc.o physconst.o
HEADER  = $(INCDIR)/body.h $(INCDIR)/neuosc.h $(INCDIR)/physconst.h $(INCDIR)/global.h

# Compilation rules
lib libneutrinocommon $(LIBNAME): $(OBJ) makefile $(HEADER)
	ar -crv $(LIBDIR)/$(LIBNAME) $(OBJ)
	ranlib $(LIBDIR)/$(LIBNAME)

%.o : %.cpp $(HEADER)
	$(GCC) -c $(CFLAGS) $(INCFLAGS) $<

.PHONY: clean
clean:
	rm -f $(LIBDIR)/$(LIBNAME) $(OBJ) 
