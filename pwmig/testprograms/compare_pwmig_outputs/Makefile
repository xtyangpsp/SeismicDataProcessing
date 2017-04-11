BIN=compare_pwmig_outputs

cxxflags=-g  -DNO_ANTELOPE
ldlibs=-lseispp -lgclgrid  -lcoords -lstock -lseispp -lperf -lm

include $(PWMIGMAKE)

OBJS=compare_pwmig_outputs.o
$(BIN) : $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
