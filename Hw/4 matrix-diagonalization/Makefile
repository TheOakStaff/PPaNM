CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all

all: out.txt

out.txt: main
	./$< > $@
	cat out.txt

main: main.c
	
