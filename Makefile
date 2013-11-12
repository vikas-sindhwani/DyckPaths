SOFTWARE_HOME=${HOME}/Software

CXX=g++

ifeq ($(shell uname -s),Darwin)
  EXTRAFLAG=-fopenmp
endif
ifeq ($(shell uname -s),Linux)
  EXTRAFLAG=-fopenmp
endif

EXTRALIB:=$(EXTRALIB) \
          -L${SOFTWARE_HOME}/lib \
	        -lboost_program_options

DBGFLAG=-g -O0
OPTIFLAG=-O3

INCLIBS = -lm $(EXTRALIB)
INCFLAGS= -I${SOFTWARE_HOME}/include
CXXFLAGS:=$(CXXFLAGS) -Wno-write-strings -Wno-format -Wno-deprecated

HEADERS = generate_Galton_Watson.hpp \
	        dyck_options.hpp \
					dyck_path.hpp

all: harness

harness-dbg: harness.cpp $(HEADERS)
	${CXX} ${CXXFLAGS} ${DBGFLAG} ${INCFLAGS} $< -o $@ ${INCLIBS}

harness: harness.cpp $(HEADERS) 
	${CXX} ${CXXFLAGS} ${OPTIFLAG} ${INCFLAGS} $< -o $@ ${INCLIBS}

harness-thr: harness.cpp $(HEADERS) 
	${CXX} ${CXXFLAGS} ${EXTRAFLAG} ${OPTIFLAG} ${INCFLAGS} $< -o $@ ${INCLIBS}

.PHONY: clean

clean:
	rm -rf *.o harness harness-dbg *.dSYM out error foo.* options.* stdout
