CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors
CFLAGS+=$(LOL)
LDFLAGS=-flto



ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4   -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=bwise

all: $(EXEC)

bwise.o: Bwise.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

bgreat: bwise.o
	$(CC) -o $@ $^ $(LDFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)

