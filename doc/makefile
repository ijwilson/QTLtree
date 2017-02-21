

all: manual.pdf

manual.pdf:  manual.tex
	R CMD pdflatex manual.tex
	R CMD pdflatex manual.tex

manual.tex: manual.Rnw
	R CMD Sweave manual.Rnw

clean:	
	rm -f *.aux manual.tex *.out *.pdf *.eps *~ *.log
