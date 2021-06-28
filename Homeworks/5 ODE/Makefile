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
	 "data.txt" using 1:2 with line title "Stepper result"\
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
	 "data2.txt" using 1:3 with line title "Infected Tc 6"\
	,"data3.txt" using 1:3 with line title "Infected Tc 3"\
	,"data4.txt" using 1:3 with line title "Infected Tc 1.5"\
	,"data5.txt" using 1:3 with line title "Infected Tc 0.75"\
	' | pyxplot

all: out.txt

out.txt: main
	./main > out.txt
	cat out.txt

main: main.o My_vec_calc.o

clean:
	$(RM) *.txt main *.png *.o
