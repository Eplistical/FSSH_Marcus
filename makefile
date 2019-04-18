CXX = g++ 
MPICXX = mpicxx
OPT = -O3
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 

all : fssh_1d afssh_1d fssh_nd_mpi afssh_1d_mpi afssh_joe_tully3_mpi afssh_landry_tully3_mpi

fssh_1d: fssh_1d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

afssh_1d: afssh_1d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_nd_mpi: fssh_nd_mpi.cpp
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

afssh_1d_mpi: afssh_1d_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

afssh_joe_tully3_mpi: afssh_joe_tully3_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

afssh_landry_tully3_mpi: afssh_landry_tully3_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

