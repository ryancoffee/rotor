CC=g++
BOOSTROOT=/opt/boost/include
CFLAGS=-Wall -I/usr/local/include -I$(BOOSTROOT) -I$(IDIR) -std=gnu++14 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local/lib -lfftw3 -lm -fopenmp

SDIR=./src
IDIR=./include
ODIR=./objects
_SRCS=Ensemble.cpp Molecules.cpp Propagators.cpp PulseFreq.cpp PulseTime.cpp ScanParams.cpp jEnsemble.cpp vEnsemble.cpp rotormain.cpp
_HEADS=Constants.hpp Ensemble.hpp Propagator.hpp PulseFreq.hpp ScanParams.hpp DataOps.hpp Molecules.hpp Propagators.hpp PulseTime.hpp jEnsemble.hpp vEnsemble.hpp rotormain.hpp

OBJECTS=$(patsubst %,$(ODIR)/%,$(_SRCS:.cpp=.o))
HEADERS=$(patsubst %,$(IDIR)/%,$(_HEADS))
SOURCES=$(patsubst %,$(SDIR)/%,$(_SRCS))
EXECUTABLE ?= rotormain

all: $(SOURCES) $(EXECUTABLE)
        
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@ 

.PHONY: clean

clean:
	rm $(ODIR)/*.o $(EXECUTABLE)

# DO NOT DELETE



