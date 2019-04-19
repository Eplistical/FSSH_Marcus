CXX = g++ 
MPICXX = mpicxx
OPT = -O3 -std=c++11
INC += 
LIBS += -lboost_program_options -lopenblas -llapack -lpthread -lgfortran 

all : afssh_1dsb_mpi

afssh_1dsb_mpi: afssh_1dsb_mpi.cpp 
	$(MPICXX) $(INC) $(OPT) $< -o $@ $(LIBS)
