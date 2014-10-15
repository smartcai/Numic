CC=gcc
INC=../
FILES=main.c ../libnumic.c
CFLAGS=-I$(INC) -lm

.PHONY: all
all:
	$(CC) $(FILES) $(CFLAGS)

.PHONY: clean
clean:
	rm *.o $(EXE)
