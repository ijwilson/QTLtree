splitfiles=	 BeagleData.o options.o RandomisedCases.o binode.o newio.o

VPATH=../sharedcode
INCLUDES=-I/home/nijw/git/QTLtree/sharedcode

#SVNDEF := -D'SVN_REV="$(shell svnversion -n .)"'
GITDEF := -D'GIT_REV="$(shell git rev-list --count --first-parent HEAD)"'

OPTFLAGS=  -O3 -DNDEBUG 
OPTFLAGS=-g  -ansi -pedantic

CXXFLAGS= $(OPTFLAGS)  $(INCLUDES) -Wall  $(GITDEF)

GSLLFLAGS=$(shell gsl-config --libs)
GSLCFLAGS= $(shell gsl-config --cflags)
LDFLAGS=  $(GSLLFLAGS) 

all:	 fastsplit QTLtree

clean:	
		rm -f QTLtree *.o 

tags:	
		etags *.cpp *.cc *.[hH] $(HOME)/utilcpp/*.cpp $(HOME)/utilcpp/*.[hH] \
		$(HOME)/GenTL/*.cpp $(HOME)/GenTL/*.[hH]

QTLtree:	 $(splitfiles) 	QTLtree.o read_csv.o  FSData.o
		$(CXX)  $(LDFLAGS) $(CXXFLAGS) $(splitfiles) FSData.o QTLtree.o  read_csv.o -o ../bin/QTLtree $(LDFLAGS)


fastsplit:	 $(splitfiles) 	fastsplit.o read_csv.o  FSData.o
		$(CXX)  $(LDFLAGS) $(CXXFLAGS) $(splitfiles) FSData.o fastsplit.o  read_csv.o -o ../bin/fastsplit  $(LDFLAGS)


RandomisedCases.o:	RandomisedCases.cpp RandomisedCases.h gsl_rand.h
fastsplit.o:	fastsplit.cc splitter.h  binode.h
FSData.o:	FSData.h read_csv.h
newio.o:	newio.cpp newio.h
QTLtree.o:		FSData.h QTLtree.cc options.h splitter.h binode.h closest.h
