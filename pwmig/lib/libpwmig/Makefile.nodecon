LIB=libpwmig.a
INCLUDE=PwmigFileHandle.h SlownessVectorMatrix.h
#MAN3=		# man pages for library routines
cxxflags=-DNO_ANTELOPE
include $(PWMIGMAKE)
#CXXFLAGS+=-D_FILE_OFFSET_BITS=64
#CXXFLAGS+=-D_LARGEFILE_SOURCE
#CXXFLAGS+=-D_LARGEFILE64_SOURCE
#  needed for debug
#CXXFLAGS+=-I/usr/OpenMotif/include
# for valgrid use only
#CXXFLAGS+=-O0

OBJS=PwmigFileHandle.o SlownessVectorMatrix.o
$(LIB) : $(OBJS)
	$(RM) $@
	$(CXXAR) $(CXXARFLAGS) $@ $(OBJS)  
	$(RANLIB) $@
