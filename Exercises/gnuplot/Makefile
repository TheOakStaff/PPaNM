CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

all: erf.pyxplot.png gam.pyxplot.png

erf.pyxplot.png: data.txt TabData.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"TabData.txt" using 1:2 with points pointtype 3 title "tabulated data"\
	' | pyxplot

gam.pyxplot.png: gamdata.txt TabData2.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set yrange [-5:5];\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "tgamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:3 with line title "homemade Gam"\
	,"TabData2.txt" using 1:2 with points pointtype 3 title "tabulated data"\
	' | pyxplot

data.txt gamdata.txt: main
	./main

main: main.o erf.o gam.o

.PHONEY: clean

clean:
	$(RM) *.o data.txt *.png main
