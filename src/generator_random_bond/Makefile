SHELL    =     /bin/bash
CFLAGS   =     -O3 -Wall -W -pedantic -Wno-deprecated
CXXFLAGS =     -O3 -Wall -W -pedantic
CXX      =     g++
CCC      =     $(CXX)
CCFLAGS  =     $(CXXFLAGS)
LIBS	 =     -lgslcblas -lgsl
PROGNAME =     ../../bin/isingGeneratorRandomBond


OBJECTS  =     main.o

all:	$(OBJECTS) executable

clean:
	rm -f $(OBJECTS) $(PROGNAME)
tidy:
	rm -f $(OBJECTS)

main.o:

executable:
	$(CXX) $(CXXFLAGS) -o $(PROGNAME) $(OBJECTS) $(LIBS)
