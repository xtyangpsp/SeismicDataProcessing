BIN=convert_ak135
cflags=-g
ldlibs=$(DBLIBS)

include $(ANTELOPEMAKE)  	
OBJS=	convert_ak135.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
