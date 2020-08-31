
.PHONY: libneutrinocommon mains

all: libneutrinocommon mains

libneutrinocommon:
	cd src; make
    
mains:
	cd mains; make

clean:
	cd src; make clean
	cd mains; make clean
	rm -f lib/*
