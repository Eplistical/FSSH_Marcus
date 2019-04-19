CXX = g++ 
MPICXX = mpicxx
OPT = -O3
INC += 
LIBS += -lboost_program_options -lopenblas -llapack -lpthread -lgfortran 

all : afssh_1dsb_mpi

afssh_1dsb_mpi: afssh_1dsb_mpi.cpp 
	$(MPICXX) $(INC) $(OPT) $< -o $@ $(LIBS)
