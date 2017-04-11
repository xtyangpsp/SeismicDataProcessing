BIN=make_earsdb
#MAN1=
ldflags=-lseispp $(DBLIBS) -lseispp -lm
include $(ANTELOPEMAKE)
OBJS=make_earsdb.o SEEDStaGeom.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
