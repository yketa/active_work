#### BUILD PARAMETERS ####

BU=build
OB=$(BU)/objects

#### COMPILATION PARAMETERS ####

CC=g++
CFLAGS=-std=gnu++11 -O3 -Wall
LDFLAGS=
MPIFLAGS=

ifneq ($(DEFINITIONS),)
	CFLAGS+=$(foreach definition, $(DEFINITIONS),-D$(definition))
endif
ifeq ($(TEST),yes)
	EXEC=$(BU)/test
	CPP=test.cpp
	LDFLAGS+=-fopenmp  # compile with openMP
	MPIFLAGS+=-fopenmp # compile with openMP
else
ifeq ($(CLONING),yes)
	EXEC=$(BU)/cloning
	CPP=cloning.cpp
	LDFLAGS+=-fopenmp  # compile with openMP
	MPIFLAGS+=-fopenmp # compile with openMP
ifeq ($(CONTROLLED_DYNAMICS),1)
	EXEC:=$(EXEC)_C1
	CFLAGS+=-DCONTROLLED_DYNAMICS=1
endif
ifeq ($(CONTROLLED_DYNAMICS),2)
	EXEC:=$(EXEC)_C2
	CFLAGS+=-DCONTROLLED_DYNAMICS=2
endif
ifeq ($(CONTROLLED_DYNAMICS),3)
	EXEC:=$(EXEC)_C3
	CFLAGS+=-DCONTROLLED_DYNAMICS=3
endif
else
ifeq ($(SIM0),yes)
	EXEC=$(BU)/simulation0
	CPP=main0.cpp
else
	EXEC=$(BU)/simulation
	CPP=main.cpp
endif
endif
endif
ifeq ($(DEBUG),yes)
	CFLAGS+=-DDEBUG
endif
ifneq ($(EXEC_NAME),)
	EXEC=$(BU)/$(EXEC_NAME)
endif
MAIN=main.cpp main0.cpp cloning.cpp test.cpp                       # files with main()
SRC=$(filter-out $(filter-out $(CPP), $(MAIN)), $(wildcard *.cpp)) # compile all files but the ones with wrong main()

ifeq ($(CELLLIST),yes)
	EXEC:=$(EXEC)_cell_list
	CFLAGS+=-DUSE_CELL_LIST
endif

ifeq ($(HEUN), yes)
	CFLAGS+=-DHEUN=true
endif

OBJ=$(addprefix $(OB)/, $(SRC:.cpp=.o))

.PHONY: all memcheck massif clean mrproper

#### COMPILATION #####

all: dir $(EXEC)

dir:
	@mkdir -p $(BU)
	@mkdir -p $(OB)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

#### DEPENDENCIES ####

$(OB)/dat.o: dat.cpp dat.hpp readwrite.hpp
	$(CC) -o $(OB)/dat.o -c dat.cpp $(CFLAGS)

$(OB)/env.o: env.cpp env.hpp
	$(CC) -o $(OB)/env.o -c env.cpp $(CFLAGS)

$(OB)/iteration.o: iteration.cpp iteration.hpp particle.hpp
	$(CC) -o $(OB)/iteration.o -c iteration.cpp $(CFLAGS)

$(OB)/maths.o: maths.cpp maths.hpp
	$(CC) -o $(OB)/maths.o -c maths.cpp $(CFLAGS)

$(OB)/particle.o: particle.cpp particle.hpp maths.hpp readwrite.hpp
	$(CC) -o $(OB)/particle.o -c particle.cpp $(CFLAGS)

$(OB)/readwrite.o: readwrite.cpp readwrite.hpp
	$(CC) -o $(OB)/readwrite.o -c readwrite.cpp $(CFLAGS)

##

$(OB)/cloning.o: cloning.cpp cloningserial.hpp env.hpp particle.hpp readwrite.hpp
	$(CC) -o $(OB)/cloning.o -c cloning.cpp $(CFLAGS) $(MPIFLAGS)

$(OB)/main.o: main.cpp env.hpp iteration.hpp particle.hpp
	$(CC) -o $(OB)/main.o -c main.cpp $(CFLAGS)

$(OB)/main0.o: main0.cpp env.hpp fire.hpp iteration.hpp maths.hpp particle.hpp
	$(CC) -o $(OB)/main0.o -c main0.cpp $(CFLAGS)

$(OB)/test.o: test.cpp
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
	rm -rf $(OB)

mrproper: clean
	rm -rf $(BU)
