all: findTauElec makeTGraphs calculateLifetimeWithScope findAllAverages findAllAveragesDigitiser decodeEvent calculateLifetimeWithDigitiser makeFFTs findFilteredAverages resizeRawAverages

PROGNAME    = findTauElec makeTGraphs findAllAverages findAllAveragesDigitiser
SOURCES     = findTauElec.cxx UsefulFunctions.cxx LifetimeConventions.cxx makeTGraphs.cxx \
RawDigitiser.cxx calculateLifetimeWithScope.cxx
INCLUDES    = LifetimeConventions.h 
OBJECTS     = $(patsubst %.cxx, %.o, $(SOURCES))
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS  := $(shell root-config --glibs)
ROOTLIBS   := $(shell root-config --nonew --libs)
CFLAGS     += $(ROOTCFLAGS)
LIBS       += $(ROOTLIBS)

LDFLAGS     = -O -fpermissive
CFLAGS  += -I$(FFTW_UTIL_INC_DIR) -fpermissive
LDFLAGS += -L$(FFTW_UTIL_INSTALL_DIR)/lib

findTauElec : findTauElec.o UsefulFunctions.o LifetimeConventions.o
	g++ -o $@ findTauElec.o UsefulFunctions.o LifetimeConventions.o $(LDFLAGS) $(LIBS) -I.

decodeEvent : decodeEvent.o UsefulFunctions.o LifetimeConventions.o
	g++ -o $@ decodeEvent.o UsefulFunctions.o LifetimeConventions.o $(LDFLAGS) $(LIBS) -I.

makeTGraphs : makeTGraphs.o RawDigitiser.o LifetimeConventions.o dict.o
	g++ -o $@ makeTGraphs.o RawDigitiser.o dict.o LifetimeConventions.o $(LDFLAGS) $(LIBS) -I.

findAllAverages : findAllAverages.o RawDigitiser.o LifetimeConventions.o UsefulFunctions.o dict.o
	g++ $(CFLAGS) -o $@ findAllAverages.o RawDigitiser.o LifetimeConventions.o UsefulFunctions.o dict.o $(LDFLAGS) $(LIBS) -lRootFftwWrapper -I.

findAllAveragesDigitiser : findAllAveragesDigitiser.o RawDigitiser.o LifetimeConventions.o UsefulFunctions.o dict.o
	g++ $(CFLAGS) -o $@ findAllAveragesDigitiser.o RawDigitiser.o LifetimeConventions.o UsefulFunctions.o dict.o $(LDFLAGS) $(LIBS) -lRootFftwWrapper -I.

calculateLifetimeWithScope : calculateLifetimeWithScope.o LifetimeConventions.o UsefulFunctions.o
	g++ -o $@ calculateLifetimeWithScope.o LifetimeConventions.o UsefulFunctions.o $(LDFLAGS) $(LIBS) -lRootFftwWrapper -I.

calculateLifetimeWithDigitiser : calculateLifetimeWithDigitiser.o LifetimeConventions.o UsefulFunctions.o 
	g++ -o $@ calculateLifetimeWithDigitiser.o LifetimeConventions.o UsefulFunctions.o $(LDFLAGS) $(LIBS) -lRootFftwWrapper -I.

findFilteredAverages : findFilteredAverages.o LifetimeConventions.o UsefulFunctions.o 
	g++ -o $@ findFilteredAverages.o LifetimeConventions.o UsefulFunctions.o $(LDFLAGS) $(LIBS) -lRootFftwWrapper -I.

resizeRawAverages : resizeRawAverages.o LifetimeConventions.o UsefulFunctions.o 
	g++ -o $@ resizeRawAverages.o LifetimeConventions.o UsefulFunctions.o $(LDFLAGS) $(LIBS) -lRootFftwWrapper -I.

makeFFTs : makeFFTs.o LifetimeConventions.o UsefulFunctions.o 
	g++ -o $@ makeFFTs.o LifetimeConventions.o UsefulFunctions.o $(LDFLAGS) $(LIBS) -lRootFftwWrapper -I.

findTauElec.o : findTauElec.cxx UsefulFunctions.cxx LifetimeConventions.h
	g++ ${CFLAGS} -c -g -o $@ $<

decodeEvent.o : decodeEvent.cxx UsefulFunctions.cxx LifetimeConventions.h
	g++ ${CFLAGS} -c -g -o $@ $<

UsefulFunctions.o : UsefulFunctions.cxx LifetimeConventions.h LifetimeConventions.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

LifetimeConventions.o : LifetimeConventions.cxx LifetimeConventions.h 
	g++ ${CFLAGS} -c -g -o $@ $<

makeTGraphs.o : makeTGraphs.cxx RawDigitiser.cxx LifetimeConventions.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

findAllAverages.o : findAllAverages.cxx LifetimeConventions.cxx UsefulFunctions.cxx RawDigitiser.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

findAllAveragesDigitiser.o : findAllAveragesDigitiser.cxx LifetimeConventions.cxx UsefulFunctions.cxx RawDigitiser.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

calculateLifetimeWithScope.o : calculateLifetimeWithScope.cxx LifetimeConventions.cxx UsefulFunctions.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

calculateLifetimeWithDigitiser.o : calculateLifetimeWithDigitiser.cxx LifetimeConventions.cxx UsefulFunctions.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

makeFFTs.o : makeFFTs.cxx LifetimeConventions.cxx UsefulFunctions.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

findFilteredAverages.o : findFilteredAverages.cxx LifetimeConventions.cxx UsefulFunctions.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

resizeRawAverages.o : resizeRawAverages.cxx LifetimeConventions.cxx UsefulFunctions.cxx
	g++ ${CFLAGS} -c -g -o $@ $<

# Generate dictionary
dict.cxx : RawDigitiser.h LinkDef.h
	rootcint -f $@ -c $(CFLAGS) -p $^

dict.o : dict.cxx RawDigitiser.cxx 
	g++ $(CFLAGS) -c -g -o $@ $<

RawDigitiser.o : RawDigitiser.cxx dict.cxx
	g++ $(CFLAGS) -c -g -o $@ $<

test:
	@echo $(ROOTCFLAGS)

clean:
	-rm -f ${PROGNAME} ${OBJECTS}
