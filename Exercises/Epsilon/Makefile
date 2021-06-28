CC = gcc # using gcc as the c-compiler
CFLAGS = -O -std=gnu11 -fwrapv # -O is the optimization, use standart gnu11
LDLIBS = -lm # lm = link math libary

out.txt: main
	./main > out.txt
	cat out.txt

main: equal.o epsilon.o
	$(CC) $^ -o $@

equal.o: equal.c
	$(CC) -c $^ -o $@ $(CFLAGS)

epsilon.o: epsilon.c
	$(CC) -c $^ -o $@ $(CFLAGS)

.PHONEY: clean

clean:
	$(RM) *.o *.txt main
