BIN=gridstacker
DIRS=
PF=gridstacker.pf

#MAN1=gridstacker.1
cxxflags=-g -DNO_ANTELOPE
ldflags=-lseispp -lgclgrid $(TRLIBS) $(DBLIBS) -lperf -lseispp -lm
ldflags=-lseispp -lgclgrid -lstock -lcoords -lperf -lseispp -lm

include $(PWMIGMAKE) 
CXXFLAGS+=-D_FILE_OFFSET_BITS=64
CXXFLAGS+=-D_LARGEFILE_SOURCE
CXXFLAGS+=-D_LARGEFILE64_SOURCE
OBJS= gridstacker.o \
  GridScratchFileHandle.o \
  GridStackPenaltyFunction.o
mdutils.cc : pwmig.h
pwmig.cc : pwmig.h
ray_transformation_operator.cc : pwmig.h
raygrid.cc : pwmig.h
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
