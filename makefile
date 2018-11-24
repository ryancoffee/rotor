CC=g++
BOOSTROOT=/opt/boost/include
SDIR=./src
IDIR=./include
ODIR=./objects
CFLAGS=-Wall -I/usr/local/include -I$(BOOSTROOT) -I$(IDIR) -std=gnu++14 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local/lib -lfftw3 -lm -fopenmp
_SRCS=Molecules.cpp jEnsemble.cpp vEnsemble.cpp PulseTime.cpp Propagators.cpp ScanParams.cpp rotormain.cpp
_HEADS=Constants.hpp DataOps.hpp Molecules.hpp PulseTime.hpp jEnsemble.hpp vEnsemble.hpp Propagators.hpp ScanParams.hpp rotormain.hpp

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

