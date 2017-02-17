SUFFIXES :=
%.c:
%.cpp:
%.o:
%.h:


CC = gcc
CPP = g++
CPPFLAGS = -std=c++11
Forall = RK.h EU.h Others.h Step.h Controllers.h ControllerADP.h Dynamical.h AlgorithmRLS.h AlgorithmADP.h AlgorithmPI.h AlgorithmVI.h MatrixCalc.h Matrix.h SquareMatrix.h SymmetricMatrix.h Diagonal.h

all: ADPsolver clean
myLua: lib clean 

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
Others.o: Others.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c Others.cpp
RK.o: RK.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c RK.cpp
EU.o: EU.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c EU.cpp

libControllerVI.o: libControllerVI.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c libControllerVI.cpp

ADPsolver: RK.o EU.o Others.o AlgorithmRLS.o AlgorithmPI.o AlgorithmVI.o Diagonal.o SymmetricMatrix.o SquareMatrix.o Matrix.o MatrixCalc.o main.o
	$(CPP) $(CPPFLAGS) -o ADPsolver RK.o EU.o Others.o AlgorithmRLS.o AlgorithmPI.o AlgorithmVI.o Diagonal.o SymmetricMatrix.o SquareMatrix.o Matrix.o MatrixCalc.o main.o

#lib: 
	#g++ -std=c++11 -shared -fPIC -o libControllerVI.so libControllerVI.cpp ControllerADP.cpp

lib: RK.o EU.o Others.o AlgorithmRLS.o AlgorithmPI.o AlgorithmVI.o Diagonal.o SymmetricMatrix.o SquareMatrix.o Matrix.o MatrixCalc.o libControllerVI.o
	$(CPP) $(CPPFLAGS) -shared -fPIC -o libControllerVI.so  RK.o EU.o Others.o AlgorithmRLS.o AlgorithmPI.o AlgorithmVI.o Diagonal.o SymmetricMatrix.o SquareMatrix.o Matrix.o MatrixCalc.o libControllerVI.o

clean:
	rm -f *.o
