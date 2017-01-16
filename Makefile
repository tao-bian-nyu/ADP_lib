SUFFIXES :=
%.c:
%.cpp:
%.o:
%.h:


CC = gcc
CPP = g++
CPPFLAGS = -std=c++11
Forall = Makefile Others.h Step.h Controllers.h ControllerADP.h Dynamical.h AlgorithmRLS.h AlgorithmADP.h AlgorithmPI.h AlgorithmVI.h MatrixCalc.h Matrix.h SquareMatrix.h SymmetricMatrix.h Diagonal.h

all: ADPsolver clean

Matrix.o: Matrix.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c Matrix.cpp
main.o: main.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c main.cpp
SquareMatrix.o: SquareMatrix.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c SquareMatrix.cpp
SymmetricMatrix.o: SymmetricMatrix.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c SymmetricMatrix.cpp
Diagonal.o: Diagonal.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c Diagonal.cpp
MatrixCalc.o: MatrixCalc.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c MatrixCalc.cpp
AlgorithmRLS.o: AlgorithmRLS.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c AlgorithmRLS.cpp
AlgorithmVI.o: AlgorithmVI.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c AlgorithmVI.cpp
AlgorithmPI.o: AlgorithmPI.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c AlgorithmPI.cpp
Dynamical.o: Dynamical.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c Dynamical.cpp
ControllerADP.o: ControllerADP.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c ControllerADP.cpp
Others.o: Others.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c Others.cpp

ADPsolver: Others.o ControllerADP.o Dynamical.o AlgorithmRLS.o AlgorithmPI.o AlgorithmVI.o Diagonal.o SymmetricMatrix.o SquareMatrix.o Matrix.o MatrixCalc.o main.o
	$(CPP) $(CPPFLAGS) -o ADPsolver Others.o ControllerADP.o Dynamical.o AlgorithmRLS.o AlgorithmPI.o AlgorithmVI.o Diagonal.o SymmetricMatrix.o SquareMatrix.o Matrix.o MatrixCalc.o main.o

clean:
	rm -f *.o
