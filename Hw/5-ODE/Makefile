CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default clean

default: Harm.pyxplot.png Pan.pyxplot.png

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

Pan.pyxplot.png: all Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set xlabel "Time(days)";\
	set ylabel "infected";\
	set yrange [0:5500000];\
	set title "Infections in DK";\
	plot \
	 "out2.txt" using 1:3 with line title "Infected Tc 5"\
	,"out3.txt" using 1:3 with line title "Infected Tc 2.5"\
	,"out4.txt" using 1:3 with line title "Infected Tc 1.25"\
	,"out5.txt" using 1:3 with line title "Infected Tc 0.625"\
	' | pyxplot

all: out.txt

out.txt: main
	./main > out.txt
	cat out.txt

main: main.c

clean:
	$(RM) *.txt main *.png
