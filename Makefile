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
	./mgsolve 8 5
	git add data
	git commit -m "new data file"

mg:
	./mgsolve 8 5

allgrid:
	./mgsolve 3 5
	./mgsolve 4 5
	./mgsolve 5 5
	./mgsolve 6 5
	./mgsolve 7 5
	./mgsolve 8 5
	
clean:
	rm -f *.o mgsolve
	rm -f ./data/Neumann/*.txt
	rm -f ./data/Dirichlet/*.txt
	
.PHONY : all clean
