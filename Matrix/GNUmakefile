#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ..
VPATH= . $(HOME) $(HOME)/lib/timer
# DIM=2
# osuffix:=$(DIM)D.o
# dsuffix:=$(DIM)D.d
# FFTWDIR = /usr/local/fftw_double
# FFTINCFLAGS = -I$(FFTWDIR)/include
system := $(shell uname)
CFLAGS = -g -Wall
BLASLIBFLAGS = -lcblas
ifeq ($(system),Darwin)
  BLASFLAGS = -framework Accelerate
endif
CFLAGS := -O3
CFLAGS += $(patsubst %, -I%, $(VPATH))
# CFLAGS += -I. -I$(HOME)/../Writers -I$(HOME)/../timer -I$(HOME)/../RectArray -I$(HOME)/../fftTools -I$(HOME)/../Hockney -DDIM=$(DIM)
LIBS:= -lconfig++ -lboost_system -lboost_filesystem -lgsl -lgslcblas -lm


#CXXNew = /Users/colella/Desktop/gcc-4.7-bin/usr/local/bin/g++ 
CXX=g++

SRCFILES:= $(wildcard $(patsubst %, %/*.cpp, $(VPATH)))
#SRCFILES += $(HOME)/../Hockney/Hockney.cpp
#SRCFILES:= $(wildcard ../RectArray/*.cpp ../../Writers/*.cpp)
OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES)) 

.PHONY: clean

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

# %$(osuffix): %.cpp GNUmakefile
# 	$(CXX) -c $(CFLAGS) $< -o $@
# 	$(CXX) -MM $(CFLAGS) $< > $*$(dsuffix)

# ../lib/libfft%D.a:$(wildcard $(HOME)/../fftTools/*.H $(HOME)/fftTools/../*.cpp)
# 	cd $(HOME)/../fftTools;make DIM=$*

LargeMatrix: GNUmakefile $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(LIBS) 

clean:
	rm  *.o *.exe *.d

-include $(OBJS:.o=.d)
