
CXXFLAGS = -fopenmp#-Wall  # put compiler settings here
# put linker settings here
CXX	= g++ -g
RM  	= rm -f 
MV	= mv
CP	= cp

all:ref shortread  longread  bin  clean   instbin 

ref:src/processRefGenome.o
	$(CXX) $(CXXFLAGS) src/processRefGenome.o   -o  processRefGenome

processRefGenome.o:src/processRefGenome.cpp
	$(CXX) $(CXXFLAGS) -c src/processRefGenome.cpp

shortread:src/processShortReads.o
	$(CXX)   src/processShortReads.o -o  processShortReads

processShortReads.o:src/processShortRead.cpp
	$(CXX) -c src/processShortRead.cpp

longread:src/processLongReads.o
	$(CXX) src/processLongReads.o  -o  processLongReads

processLongReads.o:src/processLongReads.cpp
	$(CXX) -c  src/processLongReads.cpp

bin:		
	mkdir bin

instbin:
	$(MV)  processRefGenome  processShortReads  processLongReads bin

clean:
	$(RM) src/processRefGenome.o src/processShortReads.o src/processLongReads.o  


