CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors -fopenmp
LDFLAGS=-flto -fopenmp



ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4 -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=n50 sequencesToNumbers numbersToSequences numbersFilter simulator path_counter maximal_sr crush_bulle path_to_kmer

all: $(EXEC)

n50.o: N50.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

simulator.o: simulator.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

path_counter.o: path_counter.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

crush_bulle: crush_bulle.o
	$(CC) -o $@ $^ $(LDFLAGS)

crush_bulle.o: crush_bulle.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

path_to_kmer: path_to_kmer.o
	$(CC) -o $@ $^ $(LDFLAGS)

path_to_kmer.o: path_to_kmer.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

path_counter: path_counter.o
	$(CC) -o $@ $^ $(LDFLAGS)

maximal_sr.o: Maximal_SR.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

maximal_sr: maximal_sr.o
	$(CC) -o $@ $^ $(LDFLAGS)

bwise.o: Bwise.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

bwise: bwise.o
	$(CC) -o $@ $^ $(LDFLAGS)

n50: n50.o
	$(CC) -o $@ $^ $(LDFLAGS)

simulator: simulator.o
	$(CC) -o $@ $^ $(LDFLAGS)

sequencesToNumbers.o: sequencesToNumbers.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

sequencesToNumbers: sequencesToNumbers.o
	$(CC) -o $@ $^ $(LDFLAGS)

numbersFilter.o: numbersFilter.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

improve_msr: improve_msr.o
	$(CC) -o $@ $^ $(LDFLAGS)

improve_msr.o: improve_msr.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

numbersFilter: numbersFilter.o
	$(CC) -o $@ $^ $(LDFLAGS)

numbersToSequences.o: numbersToSequences.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

numbersToSequences: numbersToSequences.o
	$(CC) -o $@ $^ $(LDFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)

