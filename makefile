# Project: mill
# Makefile created by Dev-C++ 4.9.9.2

gcc  = g++
g++   = gcc
WINDRES = windres.exe
RES  = 
OBJ  =  main.o allocation.o fileio.o diainput.o$(RES)
LINKOBJ  = main.o allocation.o fileio.o diainput.o$(RES)
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

