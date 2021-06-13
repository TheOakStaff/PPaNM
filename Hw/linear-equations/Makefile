CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt

out.txt: main
	./$< > $@

main:main.c

clean:
	$(RM) *.txt main
