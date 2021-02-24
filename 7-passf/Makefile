CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default:params main

test:
	echo $(CFLAGS)
	echo $(LDLIBS)
