BIN=getT2D
SUBDIR=/contrib
ldlibs=-lseispp -ltrvltm -lbrttutil $(TRLIBS) $(DBLIBS)  \
               -lperf -lgclgrid -lseispp

include $(ANTELOPEMAKE) 
ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKELOCAL)

CXXFLAGS += -I$(BOOSTINCLUDE)
LDFLAGS += -L$(BOOSTLIB)

CXXFLAGS += -g

OBJS= getT2D.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
