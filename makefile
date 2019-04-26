CXX = g++ 
MPICXX = mpicxx
OPT = -O3 -std=c++11
INC += 
LIBS += -lboost_program_options -lopenblas -llapack -lpthread -lgfortran 

all : afssh_1dsb_mpi afssh_2dsb_mpi afssh_2dflat_mpi afssh_2dsb_sin_mpi

afssh_1dsb_mpi: afssh_1dsb_mpi.cpp 
	$(MPICXX) $(INC) $(OPT) $< -o $@ $(LIBS)

afssh_2dsb_mpi: afssh_2dsb_mpi.cpp 
	$(MPICXX) $(INC) $(OPT) $< -o $@ $(LIBS)

afssh_2dflat_mpi: afssh_2dflat_mpi.cpp 
	$(MPICXX) $(INC) $(OPT) $< -o $@ $(LIBS)

afssh_2dsb_sin_mpi: afssh_2dsb_sin_mpi.cpp 
	$(MPICXX) $(INC) $(OPT) $< -o $@ $(LIBS)
