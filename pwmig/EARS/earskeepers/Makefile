BIN=earskeepers

cflags=-g
ldflags=$(DBLIBS) -lseispp

include $(ANTELOPEMAKE)  	
OBJS=earskeepers.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
