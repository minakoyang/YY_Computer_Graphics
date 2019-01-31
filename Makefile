## Compiler
CC=g++
## Linker
LD=$(CC)
## Flags
CPPFLAGS = -Wall -g -DLINUX
LFLAGS = -lglut -L/usr/lib -L/usr/X11R6/lib -lXi -lXmu -lGL -lGLU -lm

TARGETS = $(PROGFILES:.cpp=)

PROGFILES = \
        assn3.cpp \
        $(NULL)

targets default: $(TARGETS)

$(PROGFILES:.cpp=): assn3.o 
	$(CC) -o assn3 assn3.o ${LFLAGS}

depend :
	makedepend ${PROGFILES}
clean:
	rm -f assn3 *.o

