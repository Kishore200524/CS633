# Makefile for Assignment 1
CC = mpicc
CFLAGS = -O2 -lm

all: src

src: src.c
	$(CC) $(CFLAGS) -o src src.c

clean:
	rm -f src