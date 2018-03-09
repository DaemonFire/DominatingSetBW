#			Stage REC Hugo Hourcade 
#		Parallel programming of algorithm solving Dominating Set based on rank-width decomposition of graphs


# Documentation: http://www.gnu.org/software/make/manual/make.html
TARGET = bin/main

CC = gcc
CFLAGS = 
LDFLAGS = 

OBJECTS = $(patsubst src/%.c, obj/%.o, $(wildcard src/*.c))
HEADERS = $(wildcard include/*.h)

all: mkdir $(TARGET)

mkdir:
	mkdir -p obj

$(TARGET): $(OBJECTS)
	$(CC) $^ $(LDFLAGS) -o $@

obj/%.o: src/%.c $(HEADERS)
	$(CC) $(CFLAGS) -Iinclude -c $< -o $@

clean:
	rm -f $(TARGET)
	rm -rf obj

