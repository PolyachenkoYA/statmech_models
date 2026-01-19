CC=g++
CFLAGS=-Iinclude -O2
LDFLAGS_old=-lc -lm
LDFLAGS=-Llib -lfunc -lc -lm

sqrt: main libfunc
	$(CC) main.o $(LDFLAGS) -o sqrt

sqrt_old: main func
	$(CC) main.o func.o $(LDFLAGS) -o sqrt

main: main.c include/func.h
	$(CC) -c $(CFLAGS) main.c

libfunc: func
	ar rvs libfunc.a func.o
	mv libfunc.a lib

func: lib/func.c
	$(CC) -c $(CFLAGS) $(LDFLAGS) lib/func.c

clean:
	rm -f sqrt *.o
