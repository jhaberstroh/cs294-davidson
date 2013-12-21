#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME  := ..
VPATH= . $(HOME) $(HOME)/Matrix $(HOME)/lib/timer $(HOME)/lib/tminres-0.2/ 
#$(HOME)/lib/timer
system := $(shell uname)
CFLAGS := -g
CFLAGS += -O3
BLASLIBFLAGS = -lblas


ifeq ($(system),Darwin)
  BLASFLAGS = -framework Accelerate
  CFLAGS += -DMACINTOSH_OS
endif
# CFLAGS += -I. -I$(HOME)/../Writers -I$(HOME)/../timer -I$(HOME)/../RectArray -I$(HOME)/../fftTools -I$(HOME)/../Hockney -DDIM=$(DIM)
LIBS:= $(BLASLIBFLAGS) -llapack -llapacke -lconfig++ -lboost_system -lboost_filesystem -lgsl -lgslcblas -lm
CFLAGS += $(patsubst %, -I%, $(VPATH))


#CXXNew = /Users/colella/Desktop/gcc-4.7-bin/usr/local/bin/g++ 
CXX=g++

SRCFILES:= $(wildcard $(patsubst %, %/*.cpp, $(HOME)/Matrix $(HOME)/lib/timer $(HOME)/lib/tminres-0.2))
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

JD2: template_main.o GNUmakefile $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $< $(OBJS) $(LIBS) 

JD: JD.o GNUmakefile $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $< $(OBJS) $(LIBS) 

clean:
	rm  *.o *.exe *.d

-include $(OBJS:.o=.d)
