CC = gcc
CFLAGS = -Wall -g -lm -Wl,--stack,1677721600000
LINK = -Wl,--stack,1677721600000 #Stack Size set
EXT = exe

DEPS = solver.h
OBJ = solver.o main.o

#main: $(OBJ)
#	$(CC) $(CFLAGS) -o main solver.o main.o
main: $(OBJ)
	gcc $(CFLAGS) -o $@ $^

#main.o: main.c solver.h
#	$(CC) $(CFLAGS) -c main.c
#solver.o: solver.c solver.h
#	$(CC) $(CFLAGS) -c solver.c

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY : clean
clean :
	-rm $(OBJ) main.$(EXT)