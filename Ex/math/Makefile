CC = gcc # using gcc as the c-compiler
CFLAGS = -O -std=gnu11 # -O is the optimization, use standart gnu11
LDLIBS = -lm # lm = link math libary

default: out.txt
	cat out.txt


out.txt: math
	./math > out.txt

math: math.o
	$(CC) -o math math.o $(LDLIBS)

math.o: math.c
	$(CC) $(CFLAGS) -c math.c
