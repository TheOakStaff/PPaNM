CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

out.txt: Makefile
	./hello > out.txt

hello.o: hello.c
	$(CC) $(CFLAGS) -c hello.c -o hello.o

hello: hello.o world.o
	$(CC) $(LDFLAGS) hello.o world.o -o hello $(LDLIBS)

world.o: world.c
	$(CC) $(CFLAGS) -c $< -o $@
