BIN=geo2gclcart

SUBDIR=/contrib
ldlibs=-lm -lseispp -lgclgrid -lperf $(DBLIBS) 
ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)
CXXFLAGS += -g
CXXFLAGS += -I$(BOOSTINCLUDE)
LDFLAGS += -L$(BOOSTLIB)
OBJS=geo2gclcart.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
