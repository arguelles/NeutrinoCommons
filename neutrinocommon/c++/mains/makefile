# Compiler
GCC=g++-4.5

# Directories
LIBDIR        = ../lib
INCDIR        = ../inc
BINDIR        = ../bin

CFLAGS        = -Wall -g -O3 
INCFLAGS      = -I$(INCDIR)
LIBFLAGS      = -L$(LIBDIR) -lneutrinocommon
LDFLAGS       = -lm -lgsl -lgslcblas -lstdc++ -lrt -lpthread

all: demo

demo: $(LIBDIR)/libneutrinocommon.a demo.cpp
	$(GCC) $(INCFLAGS) -o $(BINDIR)/demo.exe demo.cpp $(LIBFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(BINDIR)/*
