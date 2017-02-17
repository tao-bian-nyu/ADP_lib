SUFFIXES :=
%.c:
%.cpp:
%.o:
%.h:


CC = gcc
CPP = g++
CPPFLAGS = -std=c++11
LUA = luajit
LIBLDIR= ../
Forall = $(LIBLDIR)RK.h $(LIBLDIR)EU.h $(LIBLDIR)Others.h $(LIBLDIR)Step.h $(LIBLDIR)Controllers.h $(LIBLDIR)ControllerADP.h $(LIBLDIR)Dynamical.h $(LIBLDIR)AlgorithmRLS.h $(LIBLDIR)AlgorithmADP.h $(LIBLDIR)AlgorithmPI.h $(LIBLDIR)AlgorithmVI.h $(LIBLDIR)MatrixCalc.h $(LIBLDIR)Matrix.h $(LIBLDIR)SquareMatrix.h $(LIBLDIR)SymmetricMatrix.h $(LIBLDIR)Diagonal.h

all: lib clean Run 

Matrix.o: $(LIBLDIR)Matrix.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)Matrix.cpp
SquareMatrix.o: $(LIBLDIR)SquareMatrix.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)SquareMatrix.cpp
SymmetricMatrix.o: $(LIBLDIR)SymmetricMatrix.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)SymmetricMatrix.cpp
Diagonal.o: $(LIBLDIR)Diagonal.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)Diagonal.cpp
MatrixCalc.o: $(LIBLDIR)MatrixCalc.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)MatrixCalc.cpp
AlgorithmRLS.o: $(LIBLDIR)AlgorithmRLS.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)AlgorithmRLS.cpp
AlgorithmVI.o: $(LIBLDIR)AlgorithmVI.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)AlgorithmVI.cpp
AlgorithmPI.o: $(LIBLDIR)AlgorithmPI.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)AlgorithmPI.cpp
Others.o: $(LIBLDIR)Others.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)Others.cpp
RK.o: $(LIBLDIR)RK.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)RK.cpp
EU.o: $(LIBLDIR)EU.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c $(LIBLDIR)EU.cpp
libControllerVI.o: libControllerVI.cpp $(Forall)
	$(CPP) $(CPPFLAGS) -c libControllerVI.cpp

lib: RK.o EU.o Others.o AlgorithmRLS.o AlgorithmPI.o AlgorithmVI.o Diagonal.o SymmetricMatrix.o SquareMatrix.o Matrix.o MatrixCalc.o libControllerVI.o
	$(CPP) $(CPPFLAGS) -shared -fPIC -o libControllerVI.so  RK.o EU.o Others.o AlgorithmRLS.o AlgorithmPI.o AlgorithmVI.o Diagonal.o SymmetricMatrix.o SquareMatrix.o Matrix.o MatrixCalc.o libControllerVI.o

clean:
	rm -f *.o
Run: 
	$(LUA) test.lua

