CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: lsf.pyxplot.png clean

lsf.pyxplot.png: all Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set xlabel "Time(days)";\
	set ylabel "ln(Activity)";\
	set title "Least Square Fit function";\
	plot \
	 "fit_data.txt" using 1:2 with line title "Our least square fit"\
	,"data.txt" using 1:2:3 with yerrorbars title "data"\
	' | pyxplot


all: out.txt

out.txt: main
	./$< > $@
	cat out.txt

main: main.c

clean:
	$(RM) main *.txt *.png
