# Compiler
CC := gcc

# Compiler options
FLAGS := -std=c99 -Wall -Werror

# Dependency
DEPS := alloc.h interp.h scheme.h integ.h stream.h

# Object files
OBJS := alloc.o interp.o scheme.o integ.o stream.o main.o

%.o: %.c $(DEPS)
	$(CC) $(FLAGS) -c $<

main.exe: $(OBJS)
	$(CC) -o $@ $^
	
clean:
	rm *.o *.exe *.xls