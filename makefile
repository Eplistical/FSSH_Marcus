CXX = g++ 
OPT = -O3
LIBS = -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 

all : fssh_1d

fssh_1d: fssh_1d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)
