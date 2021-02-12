CC = gcc
CFLAGS = -O -std=gnu11 -Wall
LDFLAGS = -I.
LDLIBS = -lm

out.txt: main
	./$< > $@

main: main.c

clean:
	$(RM) main out.txt
