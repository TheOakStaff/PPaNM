CC = gcc
CFLAGS += -pthread
LDLIBS = -lm

default:out.txt

out.txt:main
	./$< > $@
	cat out.txt

clean:
	$(RM) *.txt
