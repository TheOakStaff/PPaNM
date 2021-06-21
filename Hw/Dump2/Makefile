CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default clean

default: all

all: out.txt

out.txt: main
	./main > out.txt
	cat out.txt

main: main.o My_vec_calc.o

clean:
	$(RM) *.txt main *.png *.o
