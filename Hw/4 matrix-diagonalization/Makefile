CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: lsf.pyxplot.png clean

energy.pyxplot.png: all Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set xlabel "x";\
	set ylabel "Psi(x)";\
	set title "Quantum particle in a box";\
	plot \
	 "data2.txt" using 1:2 with line linewidth 3 title "numeric n=0"\
  ,"data2.txt" using 1:3 with line linewidth 3 title "numeric n=1"\
	,"data2.txt" using 1:4 with line linewidth 3 title "numeric n=2"\
	,"data2.txt" using 1:5 with line linewidth 1 title "analytical n=0"\
	,"data2.txt" using 1:6 with line linewidth 1 title "analytical n=1"\
	,"data2.txt" using 1:7 with line linewidth 1 title "analytical n=2"\
	' | pyxplot

all: out.txt

out.txt: main
	./$< > $@
	cat out.txt

main: main.o My_vec_calc.o

clean:
	$(RM) main *.txt *.o *.png
