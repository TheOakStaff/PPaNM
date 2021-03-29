CC = gcc # using gcc as the c-compiler
CFLAGS = -O -std=gnu11 -fwrapv # -O is the optimization, use standart gnu11
LDLIBS = -lm # lm = link math libary

default: out.txt
	cat out.txt # prints in shell

out.txt: epsilon equal # out.txt needs epsilon as argument
	./epsilon > out.txt # epsilon will print output to out.txt
	./equal >> out.txt # equal will be appended

epsilon: epsilon.o
	$(CC) -o epsilon epsilon.o $(LDLIBS) #-o epsilon is saying to build output into epsilon. epsilon.o is the source file.

epsilon.o: epsilon.c
	$(CC) $(CFLAGS) -c epsilon.c #complie epsilon.c into epsilon.o

clean:
	$(RM) *.o *.txt epsilon equal

equal: equal.c
