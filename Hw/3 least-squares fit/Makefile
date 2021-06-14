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
	 "fit_data.txt" using 1:2 with line title "My least square fit"\
	,"fit_data.txt" using 1:3 with line title "c delta- fit"\
	,"fit_data.txt" using 1:4 with line title "c delta+ fit"\
	,"data.txt" using 1:2:3 with yerrorbars title "data"\
	' | pyxplot


all: out.txt

out.txt: main
	./$< > $@
	cat out.txt

main: main.c

clean:
	$(RM) main *.txt *.png
