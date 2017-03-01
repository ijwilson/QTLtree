

all:	 zip 

zip:	update
		tar zcf QTLtree.tar.gz QTLtree

web:	zip manual
		scp index.html QTLtreeManual.pdf QTLtree.tar.gz nijw@finan.ncl.ac.uk:public_html/QTLtree/.
manual: 
		cd QTLtree/doc;make clean;make manual.pdf;
		mv QTLtree/doc/manual.pdf QTLtreeManual.pdf

update:		updateGenTL updateutil clean
		cp  ../src/FSData.cc ../src/FSData.H  QTLtree/src

updateutil:
		cp $(HOME)/utilcpp/util.H $(HOME)/utilcpp/progressBar.H $(HOME)/utilcpp/gsl_sort.H $(HOME)/utilcpp/read_csv.cc $(HOME)/utilcpp/read_csv.H $(HOME)/utilcpp/closest.H  $(HOME)/utilcpp/testStats.H $(HOME)/utilcpp/utilityfunctionals.H  $(HOME)/utilcpp/myenums.h  $(HOME)/utilcpp/newio.cpp $(HOME)/utilcpp/newio.H $(HOME)/utilcpp/gsl_rand.H $(HOME)/utilcpp/gsl_rand.cpp $(HOME)/utilcpp/options.H $(HOME)/utilcpp/options.cpp  $(HOME)/utilcpp/RandomisedCases.H $(HOME)/utilcpp/RandomisedCases.cpp QTLtree/src 

updateGenTL:
		cp  $(HOME)/GenTL/splitter.H $(HOME)/GenTL/binode.H $(HOME)/GenTL/binode.cc $(HOME)/GenTL/BeagleData.cpp $(HOME)/GenTL/BeagleData.H  QTLtree/src

clean:
		cd QTLtree/src; make clean
