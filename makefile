# Project: mill
# Makefile created by Dev-C++ 4.9.9.2

gcc  = g++
g++   = gcc
WINDRES = windres.exe
RES  = 
OBJ  =  main.o init.o fileio.o mathop.o motion.o force.o contact.o neighbour.o$(RES)
LINKOBJ  = main.o init.o fileio.o mathop.o motion.o force.o contact.o neighbour.o$(RES)
LIBS =  -lm  
INCS =  -I"include" 
CXXINCS =  -I"include" 
BIN  = test
CXXFLAGS = $(CXXINCS) -g3  
CFLAGS = -O -systype bsd43 $(INCS)  -g3
RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all-before:
all-after:
clean-custom:

all: all-before test all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(g++) $(LINKOBJ) -o "test"  $(LIBS)

main.o: main.c
	$(g++) -c main.c -o main.o $(CXXFLAGS)
init.o: init.c
	$(g++) -c init.c -o init.o $(CXXFLAGS)
fileio.o: fileio.c
	$(g++) -c fileio.c -o fileio.o $(CXXFLAGS)
mathop.o: mathop.c
	$(g++) -c mathop.c -o mathop.o $(CXXFLAGS)
contact.o: contact.c
	$(g++) -c contact.c -o contact.o $(CXXFLAGS)
force.o: force.c
	$(g++) -c force.c -o force.o $(CXXFLAGS)
motion.o: motion.c
	$(g++) -c motion.c -o motion.o $(CXXFLAGS)
neighbour.o: neighbour.c
	$(g++) -c neighbour.c -o neighbour.o $(CXXFLAGS)



