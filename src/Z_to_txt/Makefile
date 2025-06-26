SHELL    =     /bin/bash
CFLAGS   =     -m64 -O3 -Wall -W -pedantic -Wno-deprecated
CXXFLAGS =     -m64 -O3 -Wall -W -pedantic
CXX      =     g++
CCC      =     $(CXX)
CCFLAGS  =     $(CXXFLAGS)
LIBS	 =     -lgslcblas -lgsl -lgmp -lgmpxx
PROGNAME =     ../../bin/isingZToTxt


OBJECTS  =     main.o FINDmatrix.o Sample.o exp_log.o

all:	$(OBJECTS) executable

clean:
	rm -f $(OBJECTS) $(PROGNAME)
tidy:
	rm -f $(OBJECTS)

main.o:		FINDmatrix.h Sample.h exp_log.h dataType.h
FINDmatrix.o:	Sample.h exp_log.h dataType.h
Sample.o:	exp_log.h dataType.h
exp_log.o:	dataType.h
executable:
	$(CXX) $(CXXFLAGS) -o $(PROGNAME) $(OBJECTS) $(LIBS)
