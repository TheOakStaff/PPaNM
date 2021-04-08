CC = gcc # using gcc as the c-compiler
CFLAGS = -O -std=gnu11 -fwrapv # -O is the optimization, use standart gnu11
LDLIBS = -lm # lm = link math libary

out.txt: main
	./main > out.txt
	cat out.txt

main: komplex.o
	$(CC) $^ -o $@

komplex.o: komplex.c
	$(CC) -c $^ -o $@ $(CFLAGS)

.PHONEY: clean

clean:
	$(RM) *.o *.txt main
