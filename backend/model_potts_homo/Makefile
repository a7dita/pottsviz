################################################################################
#
# Makefile to compile and link C++ program, written by Kwang Il Ryom
#
# Version valid for Linux machines, based on GNU Make.
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files of .o
#
# Dependencies on included files are automatically taken care of.
#
################################################################################
# compiler
CC=g++
# flags
CFLAGS=-c -Wall

# search path for modules

INCPATH = include
INCLUDE = $(addprefix -I,$(INCPATH))


#################################################
all: mkxeq

mkxeq: rand_gen.o functions.o main.o
	$(CC) rand_gen.o functions.o main.o -o runsim
	#rm -rf *o
        
main.o: main.cpp
	$(CC) $(INCLUDE) $(CFLAGS) main.cpp

functions.o: functions.cpp
	$(CC) $(INCLUDE) $(CFLAGS) functions.cpp     

rand_gen.o: rand_gen.cpp
	$(CC) $(INCLUDE) $(CFLAGS) rand_gen.cpp 

clean:
	rm -rf *o
.PHONY: clean        
   
   
        
