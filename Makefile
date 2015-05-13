CC=g++
CFLAGS= -O3 -Wall -Winline -Wshadow -std=c++11 -fpermissive
#CFLAGS= -fpermissive
LDFLAGS=
SOURCES=mgsolve.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mgsolve
COMMON=

all: clean mgsolve

mgsolve:
	$(CC) $(CFLAGS) $(SOURCES) -o mgsolve
	
test:
	./mgsolve 10 2

clean:
	rm -f *.o mgsolve
	rm -rf ./data/arch/*.*
	mv ./data/*.* ./data/arch/

.PHONY : all clean