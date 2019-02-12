BIN=predict_firstarrivals
SUBDIR=/contrib
ldlibs=-lseispp -ltrvltm -lbrttutil $(TRLIBS) $(DBLIBS) 

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE) 
include $(ANTELOPEMAKELOCAL)

CXXFLAGS += -I$(BOOSTINCLUDE)
LDFLAGS += -L$(BOOSTLIB)

OBJS= predict_firstarrivals.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
