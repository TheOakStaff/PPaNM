CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all clean

all: int.pyxplot.png

int.pyxplot.png: my_data.txt data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Linear Interpolation";\
	plot \
	 "data.txt" using 1:2 with line title "Data"\
	,"data.txt" using 1:2 with points pointtype 3 pointsize 0.5 title "Tab value"\
	,"my_integral.txt" using 1:2 with line title "Analytical integral"\
	,"GSL_data.txt" using 1:2 with line title "GSL spline"\
	,"GSL_data.txt" using 1:3 with line title "GSL integral"\
	' | pyxplot

my_data.txt: main
	./main

main: main.c

clean:
	$(RM) *.txt main *.png
