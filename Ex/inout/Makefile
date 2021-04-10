CC = gcc # using gcc as the c-compiler
CFLAGS = -Wall -std=gnu99 # -O is the optimization, use standart gnu11
LDLIBS = -lm # lm = link math libary

.PHONEY: all clean

all: read-cmdline read-stdin read-file

clean:
	$(RM) *.o *.txt read-cmdline read-stdin read-file
