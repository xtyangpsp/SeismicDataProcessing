BIN=pwstack
PF=pwstack.pf
DATADIR=schemas/css3.0.ext
DATA=gclgdisk gclfield
INCLUDE=pwstack_reader.h

#MAN1=pwstack.1
#cxxflags=-g -DNO_ANTELOPE
cxxflags= -DNO_ANTELOPE
# Use the following for Solaris
ldlibs= -lstock -lcoords -lpwmig -lseispp -lgclgrid -lperf -lgclgrid -lseispp -lstock -lcoords
# Following is required to work with g++ 3.2.3
#ldflags=-shared -lseispp -lgclgrid -ltrvltm -lpfstream -lbrttutil $(TRLIBS) $(DBLIBS) -lperf -lgclgrid
#
# anomaly for bigred2 intel - remove for other
# -fPIE and relink with -pie
#CXX+=-fPIE

include $(PWMIGMAKE) 

pwstack_reader.cc : pwstack_reader.h
OBJS=pwstack.o \
	apsubs.o \
	aperture.o \
	pwstack_process.o \
        pwstack_reader.o \
	compute_stack_coherence.o
$(BIN) : $(OBJS)
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
