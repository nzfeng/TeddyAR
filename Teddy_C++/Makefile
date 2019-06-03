CC = g++
FLAGS = -g -std=c++11 -Wno-deprecated-declarations

INCLUDE = -isystem../ -framework OpenGL -framework GLUT
SOURCES = main.cpp helper.cpp *.cc
LIBS = 

EXENAME = main

all: $(SOURCES)
	$(CC) $(FLAGS) -o $(EXENAME) $(INCLUDE) $(LIBDIR) $(SOURCES) $(LIBS)

clean:
	rm -f *.o $(EXENAME)

.PHONY: all clean
