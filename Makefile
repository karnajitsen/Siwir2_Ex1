CC=g++
CFLAGS= -O3 -Wall -Winline -Wshadow -std=c++11 -fpermissive
#CFLAGS= -fpermissive
LDFLAGS=
SOURCES=mgsolve.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mgsolve
COMMON=

all: mgsolve

mgsolve:
	$(CC) $(CFLAGS) $(SOURCES) -o mgsolve
	
test:
	./mgsolve 3 2

clean:
	rm -f *.o mgsolve

.PHONY : all clean