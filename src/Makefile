splitfiles=	 BeagleData.o options.o RandomisedCases.o binode.o newio.o

SVNDEF := -D'SVN_REV="$(shell svnversion -n .)"'

OPTFLAGS=  -O3 -DNDEBUG 
OPTFLAGS=-g  -ansi -pedantic

CXXFLAGS= $(OPTFLAGS)  $(INCLUDES) -Wall  $(SVNDEF)

GSLLFLAGS=$(shell gsl-config --libs)
GSLCFLAGS= $(shell gsl-config --cflags)
LDFLAGS=  $(GSLLFLAGS) 

all:	 QTLtree

clean:	
		rm -f QTLtree *.o 

tags:	
		etags *.cpp *.cc *.[hH] $(HOME)/utilcpp/*.cpp $(HOME)/utilcpp/*.[hH] \
		$(HOME)/GenTL/*.cpp $(HOME)/GenTL/*.[hH]

QTLtree:	 $(splitfiles) 	QTLtree.o read_csv.o  FSData.o
		$(CXX)  $(LDFLAGS) $(CXXFLAGS) $(splitfiles) FSData.o QTLtree.o  read_csv.o -o QTLtree


fastsplit:	 $(splitfiles) 	fastsplit.o read_csv.o  FSData.o
		$(CXX)  $(LDFLAGS) $(CXXFLAGS) $(splitfiles) FSData.o fastsplit.o  read_csv.o -o fastsplit

RandomisedCases.o:	RandomisedCases.cpp RandomisedCases.H gsl_rand.H
quicksplit.o:	quicksplit.cc splitter.H  binode.H
FSData.o:	FSData.H read_csv.H
newio.o:	newio.cpp newio.H
QTLtree.o:		FSData.H QTLtree.cc options.H splitter.H binode.H closest.H
