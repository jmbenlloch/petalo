# specific names for this package

GATEMAKEDIR = $(GATE_DIR)/Makefile
GATELIBDIR = $(GATE_DIR)/lib

DICT  = petAnalysisCint
SHLIB = libpetAnalysis.so
SOURCES = $(wildcard *.cc)
HEADERS = $(filter-out petAnalysisCint.h,$(wildcard *.h))
OBJECTS = $(filter-out RunpetAnalysis.o, $(SOURCES:.cc=.o))


# include options for this package
INCFLAGS = -I.
INCFLAGS += -I$(GATE_DIR)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(GATEMAKEDIR)/Makefile.${OSNAME}

# set compiler options for ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += '-fPIC'


# call the common GNUmakefile
include $(GATEMAKEDIR)/GNUmakefile.GATE

#all: bin

### binary compilation ###

LIBS += $(shell root-config --libs)  -lCore -lRIO -lHist
LIBS += -L$(GATE_DIR)/lib/  -lGATE -lGATEIO -lGATEUtils -lpetAnalysis

LDFLAGS = -Wl,--no-as-needed

bin: lib/$(SHLIB) RunpetAnalysis.o
	@echo "<< compiling RunpetAnalysis >>"
	@$(CXX) -g $(LDFLAGS) -o ./bin/RunpetAnalysis $^ $(LIBS) 

.PHONY: clean cleanapps

clean:  cleanapps

cleanapps:
	@rm -f ./bin/RunpetAnalysis

