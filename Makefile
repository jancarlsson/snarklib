CXX = g++
CXXFLAGS = -O2 -g3 -std=c++11 -fPIC

AR = ar
RANLIB = ranlib

default :
	@echo Build options:
	@echo make autotest_bn128 LIBSNARK_PREFIX=\<path\>
	@echo make autotest_edwards LIBSNARK_PREFIX\=\<path\>
	@echo make install DIR\=\<path\>
	@echo make doc
	@echo make clean

doc :
	markdown_py -f README.html README.md -x toc -x extra --noisy

ifeq ($(DIR),)
install :
	$(error Please provide DIR, e.g. make install DIR=/usr/local/include/snarklib)
else
LIBRARY_FILES = $(shell ls *hpp *tcc | grep -v -i autotest)

# installing just copies over the template library header files
install :
	mkdir -p $(DIR)
	cp $(LIBRARY_FILES) $(DIR)
endif

clean :
	rm -f *.o autotest_bn128 autotest_edwards README.html


################################################################################
# CURVE_ALT_BN128
#

ifeq ($(LIBSNARK_PREFIX),)
autotest_bn128 :
	$(error Please provide LIBSNARK_PREFIX, e.g. make autotest_bn128 LIBSNARK_PREFIX=/usr/local)
else
CXXFLAGS_CURVE_ALT_BN128 = \
	-I$(LIBSNARK_PREFIX)/include \
	-I$(LIBSNARK_PREFIX)/include/libsnark \
	-DCURVE_ALT_BN128 -DUSE_ASM -DUSE_ADD_SPECIAL

LDFLAGS_CURVE_ALT_BN128 = \
	-L$(LIBSNARK_PREFIX)/lib \
	-Wl,-rpath $(LIBSNARK_PREFIX)/lib \
	-lgmpxx -lgmp -lprocps -lsnark

autotest_bn128 : autotest.cpp
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGS_CURVE_ALT_BN128) $< -o autotest_bn128.o
	$(CXX) -o $@ autotest_bn128.o $(LDFLAGS_CURVE_ALT_BN128)
endif


################################################################################
# CURVE_EDWARDS
#

ifeq ($(LIBSNARK_PREFIX),)
autotest_edwards :
	$(error Please provide LIBSNARK_PREFIX, e.g. make autotest_edwards LIBSNARK_PREFIX=/usr/local)
else
CXXFLAGS_CURVE_EDWARDS = \
	-I$(LIBSNARK_PREFIX)/include \
	-I$(LIBSNARK_PREFIX)/include/libsnark \
	-DCURVE_EDWARDS -DUSE_ASM -DUSE_ADD_SPECIAL

LDFLAGS_CURVE_EDWARDS = \
	-L$(LIBSNARK_PREFIX)/lib \
	-Wl,-rpath $(LIBSNARK_PREFIX)/lib \
	-lgmpxx -lgmp -lprocps -lsnark

autotest_edwards : autotest.cpp
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGS_CURVE_EDWARDS) $< -o autotest_edwards.o
	$(CXX) -o $@ autotest_edwards.o $(LDFLAGS_CURVE_EDWARDS)
endif
