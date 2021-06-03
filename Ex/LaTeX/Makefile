CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

all: exp.pyxplot.png

exp.pyxplot.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Exponential function";\
	plot \
	 "$<" using 1:2 with line linewidth 2 title "exp from math.h"\
	,"$<" using 1:3 with line linetype 2 linewidth 2 title "exp quick-and-dirty"\
	' | pyxplot

data.txt: main
	./main > data.txt
	$(MAKE) main.pdf

main.pdf: main.tex
	pdflatex $<
	pdflatex $<

clean:
	$(RM) main *.txt *.pdf *.log *.aux
