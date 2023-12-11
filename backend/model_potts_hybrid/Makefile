# Makefile to compile and link C++ program,
# written by Kwang Il Ryom.
#
# Version valid for Linux machines, based on GNU Make.
#
# Dependencies on included files are automatically taken care of, in a nasty way at the moment.
# In later version, use wild cards to include necessary header files

#################################################
# compiler
#################################################
CC=g++

#################################################
# flags
#################################################

#CFLAGS=-c -Wall -O1 -std=c++0x
#CFLAGS=-c -Wall -O1 -DDEBUG
CFLAGS=-c -Wall -fopenmp
LDFLAGS=-fopenmp

#################################################
# modules and variables
#################################################

INCPATH = include
INCLUDE = $(addprefix -I,$(INCPATH))

SOURCES=main.cpp pnet.cpp functions.cpp rand_gen.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=runsim.exe

#################################################
# instructions
#################################################

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
.cpp.o:
	$(CC) $(INCLUDE) $(CFLAGS) $< -o $@
clean:
	rm -rf *.o $(EXECUTABLE)
.PHONY: clean
