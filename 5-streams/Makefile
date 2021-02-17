CFLAGS = -o -Wall -std=gnu11
LDLIBS = -lm

.PHONEY: default
default: out.opt.txt output.txt error.txt out.file.txt out.cmdline.txt

out.cmdline.txt: cmdline input.txt
	./cmdline $(cat input.txt) > $@

out.opt.txt: cmdopt
	./cmdopt -e 0.1 -n 10 > $@
	./cmdopt -e 0.05 -n 20 >> $@
	./cmdopt -e 0.01 -n 30 >> $@

output.txt error.txt out.file.txt: input.txt stdio Makefile
	./stdio < input.txt 1> output.txt 2> error.txt


input.txt: Makefile
	cat /dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

.PHONEY: clean
clean:
	$(RM) input.txt output.txt error.txt out.file.txt stdio cmdopt cmdline
	
