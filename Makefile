CC=gcc
CFLAGS=-std=gnu99 -Wall -Wextra -O3 -ftree-vectorize -msse3
LDFLAGS=
SOURCES=*.c
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) -DHELMHOLTZ $(CFLAGS) $(LDFLAGS) $(OBJECTS) -o $@ -lm

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) *.o