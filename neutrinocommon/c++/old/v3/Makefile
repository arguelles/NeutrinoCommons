# Compiler
GCC=g++-4.5
#GCC=gcc-4.5

# Directories
MLINKDIR      = /usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit
CADDSDIR      = $(MLINKDIR)/Linux-x86-64/CompilerAdditions/

CFLAGS        = -Wall -g -O3 -lm -lgsl -lgslcblas
INCFLAGS      = -I$(CADDSDIR)
LIBFLAGS      = -L$(CADDSDIR)
LDFLAGS       =

# Project files
BIN    = neutrinocommon
OBJ    = body.o neuosc.o physconst.o nc_MLinterface.o
TM     = neutrinocommon.tm
TM_C   = neutrinocommon-tm.c
HEADER = body.h neuosc.h physconst.h nc_MLinterface.h

# Compilation rules
$(BIN): $(OBJ)
	mcc $(CFLAGS) $(INCFLAGS) $(OBJ) $(LIBFLAGS) $(LDFLAGS) -o $(BIN)

%.o : %.cpp $(HEADER)
	$(GCC) -c $(CFLAGS) $(INCFLAGS) $<

$(TM_C) : $(TM)
	$(CADDSDIR)/mprep -prototypes $? -o $@

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ) $(TM_C)

