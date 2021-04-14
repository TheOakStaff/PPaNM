CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

erf.pyxplot.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Errro function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	' | pyxplot

data.txt: main
	./main > data.txt

main: main.o erf.o

.PHONEY: clean

clean:
	$(RM) *.o *.txt .png main
