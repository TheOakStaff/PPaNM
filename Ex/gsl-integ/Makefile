CC = gcc
# This is neede to eneable gsl fuinctions
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

Default: out.txt

out.txt: main
	./main > out.txt
	cat $@

main: main.o

main.o: main.c
