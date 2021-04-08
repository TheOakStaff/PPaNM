#CC = gcc # using gcc as the c-compiler
CFLAGS = -Wall -std=gnu99 # -O is the optimization, use standart gnu11
LDLIBS = -lm # lm = link math libary

.PHONEY: all clean

all: out.txt
	cat $<

out.txt: main
	./$< > $@

main: main.o komplex.o

main.o komplex.o:komplex.h

clean:
	$(RM) *.o *.txt main
