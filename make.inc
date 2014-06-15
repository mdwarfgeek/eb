### Constants: edit to suit your system ###

# C compiler
#CC=gcc

# Optimization flags

# DEBUG:
#OPT=-g

# OPT:
OPT=-g -O3 -ffast-math

# Compiler flags
CFLAGS=-std=gnu99 -Wall $(OPT) -fPIC

### End constants section ###

# Rules for building

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

