CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt

out.txt: main
	./$< > $@
	cat out.txt

main: main.o My_vec_calc.o SVD.o

clean:
	$(RM) main out.txt TimeTest.txt *.o
