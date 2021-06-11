CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all clean

all: int.pyxplot.png

int.pyxplot.png: my_data.txt out.txt data.txt TabData.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Interpolation of x = 2.5";\
	plot \
	 "data.txt" using 1:2 with line title "Data"\
	,"my_data.txt" using 1:2 with points pointtype 3 title "Calculated"\
	,"TabData.txt" using 1:2 with points pointtype 3 pointsize 0.5 title "Tab value"\
	' | pyxplot

my_data.txt: main

out.txt: main
	./main
	cat out.txt

main: main.c

clean:
	$(RM) out.txt my_data.txt main
