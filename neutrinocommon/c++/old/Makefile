# Compiler
GCC=g++-4.5

# Directories
CADDSDIR      = ./

CFLAGS        = -Wall -g -O3 
INCFLAGS      = -I$(CADDSDIR)
LIBFLAGS      = -L$(CADDSDIR)
LDFLAGS       = -lm -lgsl -lgslcblas -lstdc++ -lrt -lpthread

# Project files
BIN    = main.exe
OBJ    = body.o neuosc.o physconst.o
HEADER = body.h neuosc.h physconst.h global.h
MAIN   = main.cpp

# Compilation rules
$(BIN): $(OBJ)
	$(GCC) $(MAIN) $(CFLAGS) $(INCFLAGS) $(OBJ) $(LIBFLAGS) $(LDFLAGS) -o $(BIN)	

%.o : %.cpp $(HEADER)
	$(GCC) -c $(CFLAGS) $(INCFLAGS) $<

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ) 
