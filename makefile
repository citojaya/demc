# Project: mill
# Makefile created by Dev-C++ 4.9.9.2

gcc  = g++
g++   = gcc
WINDRES = windres.exe
RES  = 
OBJ  =  main.o allocation.o fileio.o diainput.o math_operation.o motion.o$(RES)
LINKOBJ  = main.o allocation.o fileio.o diainput.o math_operation.o motion.o $(RES)
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
allocation.o: allocation.c
	$(g++) -c allocation.c -o allocation.o $(CXXFLAGS)
fileio.o: fileio.c
	$(g++) -c fileio.c -o fileio.o $(CXXFLAGS)
diainput.o: diainput.c
	$(g++) -c diainput.c -o diainput.o $(CXXFLAGS)
math_operation.o: math_operation.c
	$(g++) -c math_operation.c -o math_operation.o $(CXXFLAGS)
motion.o: motion.c
	$(g++) -c motion.c -o motion.o $(CXXFLAGS)

