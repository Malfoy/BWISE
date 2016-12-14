CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors
CFLAGS+=$(LOL)
LDFLAGS=-flto



ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4 -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=bwise n50

all: $(EXEC)

n50.o: N50.cpp
	 $(CC) -o $@ -c $< $(CFLAGS)

bwise.o: Bwise.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

bwise: bwise.o
	$(CC) -o $@ $^ $(LDFLAGS)

n50: n50.o
	$(CC) -o $@ $^ $(LDFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)

