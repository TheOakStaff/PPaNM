CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: Harm.pyxplot.png clean

Harm.pyxplot.png: all Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set xlabel "Time";\
	set ylabel "y(t+h)";\
	set title "ODE 12stepper";\
	plot \
	 "out.txt" using 1:2 with line title "Stepper result"\
	' | pyxplot

all: out.txt

out.txt: main
	./main > out.txt
	cat out.txt

main: main.c

clean:
	$(RM) *.txt main *.png
