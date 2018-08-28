CC=g++
CFLAGS=-Wall -I/usr/local/include -I$(HOME)/computing/boost -I./include -I$(HOME)/cpp/include -std=gnu++11 -DHAVE_INLINE -O3 -fopenmp
LDFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lm -fopenmp

EXECUTABLE=run_rotor
SDIR=./src
IDIR=./include
ODIR=./objdir

vpath %.hpp $(IDIR)
vpath %.o $(ODIR)
vpath %.cpp $(SDIR)

DEPS := $(wildcard $(IDIR)/*.hpp)
SRCS := $(wildcard $(SDIR)/*.cpp)
_OBJS := $(patsubst %.cpp,%.o,$(wildcard $(SDIR)/*.cpp))
OBJS := $(patsubst $(SDIR)/%,$(ODIR)/%,$(_OBJS))

print-%  : ; @echo $* = $($*)

$(ODIR)/%.o : $(SDIR)/%.cpp $(DEPS)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) -c $<

$(OBJS): | $(ODIR)

$(ODIR):
	mkdir $(ODIR)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXECUTABLE)


# DO NOT DELETE

