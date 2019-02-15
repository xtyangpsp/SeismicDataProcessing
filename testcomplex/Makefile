BIN=testcomplex
#DATADIR=schemas/css3.0.ext
#DATA=decon
SUBDIR=/contrib

#PF=trace_decon.pf

ldlibs=-L$(XMOTIFLIB) -L$(BOOSTLIB) \
	-lgclgrid -lseispp -lperf -L/usr/local/Cellar/gsl/2.5/lib -lgsl -lgslcblas -lseispp $(TRLIBS) $(DBLIBS)  \
       $(F77LIBS) $(X11LIBS)

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)
CXX=g++-8
CXXFLAGS += -I$(BOOSTINCLUDE) -I/usr/local/Cellar/gsl/2.5/include -fopenmp -O3
#LDLIBS += -L/usr/local/include
LDFLAGS += -Wl,-rpath,/usr/local/Cellar/gsl/2.5/lib #/usr/local/lib

OBJS= trace_decon.o ComplexArray.o dpss.o DeconOperator.o 
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
