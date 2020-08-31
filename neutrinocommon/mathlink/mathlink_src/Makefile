# Compiler
GCC=g++-4.5

# Directories
NCMAINDIR     = /home/carlos/workspace/NeutrinoCommons/neutrinocommon
MLINKDIR      = /usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit
SYS           = Linux-x86-64
CADDSDIR      = $(MLINKDIR)/$(SYS)/CompilerAdditions

NCINCDIR	  = $(NCMAINDIR)/c++/inc
NCLIBDIR	  = $(NCMAINDIR)/c++/lib

# ML libraty. Set to ML32i3 if using a 32-bit system.
MLLIB         = ML64i3 
CFLAGS        = -Wall -g -O3 
INCFLAGS      = -I$(CADDSDIR) -I$(NCINCDIR)
LIBFLAGS      = -L$(CADDSDIR) -L$(NCLIBDIR)
LDFLAGS       = -lm -lgsl -lgslcblas -lML64i3 -lstdc++ -lrt -lpthread -lneutrinocommon

# Project files
BIN    = neutrinocommon.exe
TM     = neutrinocommon.tm
TM_C   = neutrinocommon-tm.c

# Compilation rules
$(BIN): $(TM_C) $(NCLIBDIR)/libneutrinocommon.a
	$(GCC) $(TM_C) $(CFLAGS) $(INCFLAGS) $(OBJ) $(LIBFLAGS) $(LDFLAGS) -o $(BIN)    

%.o : %.cpp $(HEADER)
	$(GCC) -c $(CFLAGS) $(INCFLAGS) $<

$(TM_C) : $(TM)
	$(CADDSDIR)/mprep -prototypes $? -o $@  

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ) $(TM_C) neutrinocommon.tm.o
