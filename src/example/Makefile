CC=gcc
INC=../
FILES=test.c ../libmatrix.c
CFLAGS=-I$(INC) -lm

.PHONY: all
all:
	$(CC) $(FILES) $(CFLAGS)

.PHONY: clean
clean:
	rm *.o $(EXE)
