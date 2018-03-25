CPPFLAGS=-Wall -Wextra -std=gnu11
CFLAGS+=-O3 -ggdb #-fopenmp

EXECUTABLES=llmco fasta2co fastq2co readco

.DELETE_ON_ERROR:
.PHONY: clean all format dist

all: $(EXECUTABLES)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $^

co2dist2: co2dist2.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^

co2dist_p: co2dist_p.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^

fasta2co: fasta2co_v18.3.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^

fastq2co: fastq2co.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^

readco: readco.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^


clean:
	rm -f $(EXECUTABLES) *.o
