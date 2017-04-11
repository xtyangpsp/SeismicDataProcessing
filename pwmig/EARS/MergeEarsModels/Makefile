BIN=MergeEarsModels

#MAN1=

cflags= -g
ldlibs= $(DBLIBS)

include $(ANTELOPEMAKE)  	

OBJS=MergeEarsModels.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
