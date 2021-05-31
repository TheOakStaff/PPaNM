CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp lgomp
LDLIBS = -lm

out.txt: main
	./$< > $@
	cat out.txt
	
clean:
	$(RM) main out.txt
