BIN=trace_decon
DATADIR=schemas/css3.0.ext
DATA=decon wfprocess sclink evlink
SUBDIR=/contrib

PF=trace_decon.pf

ldlibs= -lseispp -L$(BOOSTLIB) \
	-lperf -L/usr/local/Cellar/gsl/2.5/lib -lgsl \
	-lgslcblas $(TRLIBS) $(DBLIBS) $(F77LIBS) \
	-L/usr/local/Cellar/libomp/7.0.0/lib -lomp 
	 
ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)
include $(ANTELOPEMAKELOCAL)
#CXX=g++-8
CXXFLAGS += -I/usr/local/Cellar/libomp/7.0.0/include \
	-I/usr/local/Cellar/gsl/2.5/include -I$(BOOSTINCLUDE) -Xpreprocessor -fopenmp -O3
#LDLIBS += -L/usr/local/Cellar/gcc/8.2.0/lib
#LDFLAGS += -Wl,-rpath,/usr/local/Cellar/gsl/2.5/lib #/usr/local/lib

OBJS= trace_decon.o ComplexArray.o DeconOperator.o dpss.o dmatrix.o 
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
