#VPATH=$(HOME)/utilcpp:$(HOME)/GenTL

vpath % $(HOME)/utilcpp:$(HOME)/GenTL

splitfiles=	 BeagleData.o options.o RandomisedCases.o binode.o newio.o

#STL=stlport
ifeq	($(STL),stlport)
	STLPORTINC=-I/usr/include/stlport
	STLLIB= -lstlport
else
	STLPORTINC=
	STLLIB=
endif

INCLUDES=  -I$(HOME)/GenTL -I$(HOME)/utilcpp 

SVNDEF := -D'SVN_REV="$(shell svnversion -n .)"'

HOSTNAME := "$(shell hostname)"
ifeq ($(HOSTNAME),"masternode")
	EXFLAGS= -m32
endif

OPTFLAGS= $(EXFLAGS) -O3 -DNDEBUG  #  -fomit-frame-pointer -funroll-loops 
#OPTFLAGS= -pg -O4 -DNDEBUG  -funroll-loops
#OPTFLAGS= -O2  -Wall  -funroll-loops -march=pentium4 -mfpmath=sse -DNDEBUG   #-msse2 -mmmx 
#OPTFLAGS=-g -DTNT_BOUNDS_CHECK -ansi -pedantic
#CXXFLAGS= $(OPTFLAGS)  $(INCLUDES) -Wall -m32 $(SVNDEF)
CXXFLAGS= $(OPTFLAGS)  $(INCLUDES) -Wall  $(SVNDEF)
#CXX=g++4
AR=ar

GSLLFLAGS=$(shell gsl-config --libs)
GSLCFLAGS= $(shell gsl-config --cflags)
#LDFLAGS= $(GSLLFLAGS) -m32 $(STLLIB)
LDFLAGS= $(EXFLAGS) $(GSLLFLAGS)  $(STLLIB)
#CXX=mpic++

all:	 fastsplit fastsplitqtl QTLTree ExtractVars AlleleCounts

clean:	
	rm -f mutsim *.o 

tags:	
		etags *.cpp *.cc *.[hH] $(HOME)/utilcpp/*.cpp $(HOME)/utilcpp/*.[hH] \
		$(HOME)/GenTL/*.cpp $(HOME)/GenTL/*.[hH]

ExtractVars:  ExtractVars.o FSData.o  read_csv.o options.o newio.o 
		$(CXX)  $(LDFLAGS) $(CXXFLAGS)   ExtractVars.o FSData.o  newio.o  options.o read_csv.o -o bin/ExtractVars $(LDFLAGS)

AlleleCounts:  AlleleCounts.o FSData.o  read_csv.o options.o newio.o 
		$(CXX)  $(LDFLAGS) $(CXXFLAGS)  AlleleCounts.o FSData.o  newio.o  options.o read_csv.o -o bin/AlleleCounts $(LDFLAGS)


QTLTree:	 $(splitfiles) 	QTLTree.o read_csv.o  FSData.o
		$(CXX)  $(LDFLAGS) $(CXXFLAGS) $(splitfiles) FSData.o QTLTree.o  read_csv.o -o bin/QTLTree $(LDFLAGS)


fastsplit:	fastsplit.o $(splitfiles) 
		$(CXX)   $(CXXFLAGS) $(splitfiles) fastsplit.o  -o bin/fastsplit 	$(LDFLAGS)

fastsplitqtl:	fastsplitqtl.o $(splitfiles) 
		$(CXX)   $(CXXFLAGS) $(splitfiles) fastsplitqtl.o  -o bin/fastsplitqtl $(LDFLAGS)

fastsplit.o:	fastsplit.cc options.H splitter.H binode.H closest.H
fastsplitqtl.o:	fastsplitqtl.cc options.H splitter.H binode.H closest.H


RandomisedCases.o:	RandomisedCases.cpp RandomisedCases.H gsl_rand.H
quicksplit.o:	quicksplit.cc splitter.H  binode.H
FSData.o:	FSData.H read_csv.H
newio.o:	newio.cpp newio.H
QTLTree.o:		FSData.H QTLTree.cc options.H splitter.H binode.H closest.H
