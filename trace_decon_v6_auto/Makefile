BIN=trace_decon
DATADIR=schemas/css3.0.ext
DATA=decon
SUBDIR=/contrib

PF=trace_decon.pf

ldlibs=-L$(XMOTIFLIB) -L$(BOOSTLIB) -lseisw -lXm -lXt \
	-lgclgrid -lseispp -lperf -L/usr/local/Cellar/gsl/2.5/lib -lgsl -lgslcblas -lseispp $(TRLIBS) $(DBLIBS)  \
       $(F77LIBS) $(X11LIBS)

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)
CXX=c++-6
CXXFLAGS += -I$(BOOSTINCLUDE) -I/usr/local/Cellar/gsl/2.5/include -fopenmp -O3
LDLIBS += -L~/SOFT/usr/local/include
LDFLAGS += -Wl,-rpath,~/SOFT/usr/local/lib

OBJS= trace_decon.o ComplexArray.o DeconOperator.o dpss.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
