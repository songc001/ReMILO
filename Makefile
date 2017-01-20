
CXXFLAGS = -fopenmp#-Wall  # put compiler settings here
# put linker settings here
CXX	= g++ -g
RM  	= rm -f 
MV	= mv
CP	= cp

all:ref shortread  longread  bin  clean   instbin 

ref:src/ref.o
	$(CXX) $(CXXFLAGS) src/ref.o   -o ref

ref.o:src/ref.cpp
	$(CXX) $(CXXFLAGS) -c src/ref.cpp

shortread:src/shortread.o
	$(CXX)   src/shortread.o -o shortread

shortread.o:src/shortread.cpp
	$(CXX) -c src/shortread.cpp

longread:src/longread.o
	$(CXX) src/longread.o  -o  longread

longread.o:src/longread.cpp
	$(CXX) -c  src/longread.cpp

bin:		
	mkdir bin

instbin:
	$(MV)  ref shortread  longread bin

clean:
	$(RM) src/ref.o src/shortread.o src/longread.o  


