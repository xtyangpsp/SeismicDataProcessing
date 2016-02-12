BIN=RFrotate
MAN1=RFrotate.1
DATADIR=schemas/css3.0.ext
DATA=tmatrix
CXXFLAGS=-g
PF=RFrotate.pf
SUBDIR= /contrib

ldlibs=-lseispp -lgclgrid -lperf -ltrvltm $(DBLIBS) $(TRLIBS) -lperf -lm

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)
include $(ANTELOPEMAKELOCAL)

CXXFLAGS += -I$(BOOSTINCLUDE)

OBJS=RFrotate.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
