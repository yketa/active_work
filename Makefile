#### BUILD PARAMETERS ####

BU=build
OB=$(BU)/objects

#### COMPILATION PARAMETERS ####

CC=g++
CFLAGS=-std=gnu++11
LDFLAGS=

ifeq ($(TEST),yes)
	EXEC=$(BU)/test
	SRC=$(filter-out main.cpp, $(wildcard *.cpp))
else
	EXEC=$(BU)/simulation
	SRC=$(filter-out test.cpp, $(wildcard *.cpp))
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

$(OB)/env.o: env.cpp env.h
	$(CC) -o $(OB)/env.o -c env.cpp $(CFLAGS)

$(OB)/iteration.o: iteration.cpp iteration.h particle.h
	$(CC) -o $(OB)/iteration.o -c iteration.cpp $(CFLAGS)

$(OB)/maths.o: maths.cpp maths.h
	$(CC) -o $(OB)/maths.o -c maths.cpp $(CFLAGS)

$(OB)/particle.o: particle.cpp particle.h maths.h
	$(CC) -o $(OB)/particle.o -c particle.cpp $(CFLAGS)

$(OB)/read.o: read.cpp read.h
	$(CC) -o $(OB)/read.o -c read.cpp $(CFLAGS)

$(OB)/work.o: work.cpp work.h read.h
	$(CC) -o $(OB)/work.o -c work.cpp $(CFLAGS)

##

$(OB)/main.o: main.cpp env.h iteration.h maths.h particle.h read.h work.h
	$(CC) -o $(OB)/main.o -c main.cpp $(CFLAGS)

$(OB)/test.o: test.cpp env.h iteration.h maths.h particle.h read.h work.h
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

