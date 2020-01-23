CC = gcc
CFLAGS = -Wall -g  #-Wl,--stack,1677721600000
STACK = #-Wl,--stack,1677721600000
LINK = -lm #Stack Size set

EXT = 

DEPS = 2Dsolver.h
OBJ = 2Dsolver.o main.o

#main: $(OBJ)
#	$(CC) $(CFLAGS) -o main 2Dsolver.o main.o
main: $(OBJ)
	gcc $(CFLAGS) $(STACK) -o $@ $^ $(LINK)

#main.o: main.c 2Dsolver.h
#	$(CC) $(CFLAGS) -c main.c
#2Dsolver.o: sDolver.c 2Dsolver.h
#	$(CC) $(CFLAGS) -c 2Dsolver.c

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY : clean
clean :
	-rm $(OBJ) main