#### BUILD PARAMETERS ####

BU=build
OB=$(BU)/objects

#### COMPILATION PARAMETERS ####

CC=g++
CFLAGS=-std=gnu++11 -O3 -Wall
LDFLAGS=

ifeq ($(TEST),yes)
	EXEC=$(BU)/test
	SRC=$(filter-out main.cpp cloning.cpp, $(wildcard *.cpp))
else
ifeq ($(CLONING),yes)
	EXEC=$(BU)/cloning
	SRC=$(filter-out main.cpp test.cpp, $(wildcard *.cpp))
else
	EXEC=$(BU)/simulation
	SRC=$(filter-out test.cpp cloning.cpp, $(wildcard *.cpp))
endif
endif

ifeq ($(CELLLIST),yes)
	CFLAGS+= -DUSE_CELL_LIST
endif

OBJ=$(addprefix $(OB)/, $(SRC:.cpp=.o))

.PHONY: all clean mrproper

#### COMPILATION #####

all: dir $(EXEC)

dir:
	@mkdir -p $(BU)
	@mkdir -p $(OB)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

#### DEPENDENCIES ####

$(OB)/env.o: env.cpp env.hpp
	$(CC) -o $(OB)/env.o -c env.cpp $(CFLAGS)

$(OB)/iteration.o: iteration.cpp iteration.hpp particle.hpp
	$(CC) -o $(OB)/iteration.o -c iteration.cpp $(CFLAGS)

$(OB)/maths.o: maths.cpp maths.hpp
	$(CC) -o $(OB)/maths.o -c maths.cpp $(CFLAGS)

$(OB)/param.o: param.cpp param.hpp
	$(CC) -o $(OB)/param.o -c param.cpp $(CFLAGS)

$(OB)/particle.o: particle.cpp particle.hpp maths.hpp
	$(CC) -o $(OB)/particle.o -c particle.cpp $(CFLAGS)

$(OB)/read.o: read.cpp read.hpp
	$(CC) -o $(OB)/read.o -c read.cpp $(CFLAGS)

$(OB)/write.o: write.cpp write.hpp
	$(CC) -o $(OB)/write.o -c write.cpp $(CFLAGS)

##

$(OB)/cloning.o: cloning.cpp cloningserial.hpp env.hpp param.hpp particle.hpp
	$(CC) -o $(OB)/cloning.o -c cloning.cpp $(CFLAGS)

$(OB)/main.o: main.cpp env.hpp iteration.hpp maths.hpp param.hpp particle.hpp read.hpp
	$(CC) -o $(OB)/main.o -c main.cpp $(CFLAGS)

$(OB)/test.o: test.cpp env.hpp iteration.hpp maths.hpp param.hpp particle.hpp read.hpp
	$(CC) -o $(OB)/test.o -c test.cpp $(CFLAGS)

#### VALGRIND ####

memcheck: dir $(OBJ)
	$(CC) -g -o $(EXEC) $(OBJ) $(LDFLAGS)
	valgrind --leak-check=yes --track-origins=yes --log-file=$(BU)/memcheck.output $(EXEC)

massif: dir $(OBJ)
	$(CC) -g -o $(EXEC) $(OBJ) $(LDFLAGS)
	valgrind --tool=massif --massif-out-file=$(BU)/massif.out $(EXEC)
	ms_print $(BU)/massif.out > $(BU)/massif.output

#### CLEAN ####

clean:
	@rm -rf $(OB)

mrproper: clean
	@rm -rf $(BU)
