CC = gcc
CFLAGS = -Wall -O2 -Iinclude
LIBS = -lgsl -lgslcblas -lm

SRC = src/main.c src/solucionar.c src/potential.c src/matrix.c
OBJ = $(SRC:.c=.o)

all: schrodinger

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

schrodinger: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f src/*.o schrodinger -rf data
