CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default clean

default: Gauss.pyxplot.png

Gauss.pyxplot.png: all Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gauss";\
	plot \
	 "data.txt" using 1:2 with line title "Minimizer results"\
	,"data.txt" using 1:2 with points pointtype 3 pointsize 0.25 title "points"\
	,"FindPoint.txt" using 1:2 with points pointtype 2 pointsize 0.5 title "Random point"\
	' | pyxplot

#	,"Tabdata.txt" using 1:2 with line title "Gauss function"
#	,"Tabdata.txt" using 1:2 with points pointtype 2 pointsize 0.5 title "Gauss tab values"



all: out.txt

out.txt: main
	./main > out.txt
	cat out.txt

main: main.o vec_calc.o minimization.o

clean:
	$(RM) *.txt main *.png *.o
