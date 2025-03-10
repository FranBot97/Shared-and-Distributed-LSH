CXX = g++
MPI_CXX = mpicxx
CXXFLAGS = -std=c++20 -Wall -I. -Idependencies -O3 -DNDEBUG -ffast-math

# File sorgente principali
SEQ_SRC = sequential.cpp
FF_SRC = fastflow.cpp
MPI1_SRC = mpi1.cpp
MPI2_SRC = mpi2.cpp

# File oggetto delle dipendenze
DEPENDENCY_OBJECTS = dependencies/frechet_distance.o \
                     dependencies/frechet_distance2.o \
                     dependencies/geometry_basics.o

# Nome dei target
SEQ_TARGET = sequential
FF_TARGET = fastflow
MPI1_TARGET = mpi1
MPI2_TARGET = mpi2

# Regole per compilare 
sequential: $(SEQ_TARGET)
fastflow: $(FF_TARGET)
mpi1: $(MPI1_TARGET)
mpi2: $(MPI2_TARGET)

all: $(SEQ_TARGET) $(FF_TARGET) $(MPI1_TARGET) $(MPI2_TARGET) 

# Regola per compilare l'eseguibile sequential
$(SEQ_TARGET): $(SEQ_SRC) $(DEPENDENCY_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(SEQ_TARGET) $(SEQ_SRC) $(DEPENDENCY_OBJECTS)

# Regola per compilare l'eseguibile fastflow
$(FF_TARGET): $(FF_SRC) $(DEPENDENCY_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(FF_TARGET) $(FF_SRC) $(DEPENDENCY_OBJECTS)

# Regola per compilare l'eseguibile mpi1
$(MPI1_TARGET): $(MPI1_SRC) $(DEPENDENCY_OBJECTS)
	$(MPI_CXX) $(CXXFLAGS) -o $(MPI1_TARGET) $(MPI1_SRC) $(DEPENDENCY_OBJECTS)

# Regola per compilare l'eseguibile mpi2
$(MPI2_TARGET): $(MPI2_SRC) $(DEPENDENCY_OBJECTS)
	$(MPI_CXX) $(CXXFLAGS) -fopenmp -o $(MPI2_TARGET) $(MPI2_SRC) $(DEPENDENCY_OBJECTS)

dependencies/%.o: dependencies/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Regola per pulire tutti i file generati
clean:
	rm -f $(SEQ_TARGET) $(FF_TARGET) $(MPI1_TARGET)  $(MPI2_TARGET)  dependencies/*.o