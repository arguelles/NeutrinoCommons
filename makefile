TOP=$(shell pwd)

NUSIGMA=neutrinocommon/linked/nusigma
HEPMC_src=neutrinocommon/linked/hepmc/src 
HEPMC_bld=neutrinocommon/linked/hepmc/build
TAUOLA=neutrinocommon/linked/tauola/src

all:
	cd $(NUSIGMA); make;
	cd $(TOP)/$(HEPMC_bld); cmake ../src -Dmomentum=GEV -Dlength=CM; make install;
	cd $(TAUOLA); ./configure --with-hepmc=$(TOP)/$(HEPMC_bld); make install;
	cd neutrinocommon/shareobjects/cython_codes/taudecay; python setup.py build_ext --inplace; cp taudecay.so ../../
	python setup.py install --user;

clean:
	cd $(NUSIGMA); make clean;
	cd $(TAUOLA); make clean;
	rm -rf $(HEPMC_bld)/*;
	rm -rf build
	cd neutrinocommon/shareobjects/cython_codes/taudecay; rm taudecay.so taudecay.cpp
