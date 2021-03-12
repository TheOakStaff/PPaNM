CC = gcc # using gcc as the c-compiler
CFLAGS = -O -std=gnu11 # -O is the optimization, use standart gnu11
LDLIBS = -lm # lm = link math libary

default: out.txt
	cat out.txt # prints in shell

out.txt: hello # out.txt needs hello as argument
	./hello > out.txt # hello will print output to out.txt

hello: hello.o
	$(CC) -o hello hello.o $(LDLIBS) #-o hello is saying to build output into hello. hello.o is the source file.

hello.o: hello.c
	$(CC) $(CFLAGS) -c hello.c #complie hello.c into hello.o

.PHONEY: test
test: # Debugging
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)
