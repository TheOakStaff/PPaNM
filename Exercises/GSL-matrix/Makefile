CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

default: out.txt

out.txt: main
	./$< > $@
	cat out.txt

main: main.c

clean:
	$(RM) main *.txt
