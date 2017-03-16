ADPlib - A C++ library supporting the implementation of various ADP algorithms for users and developers from different background.
------------------------------------------------

This library contains numerical adaptive dynamic programming (ADP) algorithms to solve linear quadratic regulator (LQR) problems and algebraic Riccati equations (ARE). ADP aims at finding a stabilizing optimal control policy for dynamical systems via online learning. The aim of this project is to provide a complete implementation for online and offline ADP solvers that can be easily applied in the controller design for dynamical systems.

Learn more about LQR at: https://en.wikipedia.org/wiki/Linear–quadratic_regulator

Learn more about ARE at: https://en.wikipedia.org/wiki/Algebraic_Riccati_equation

Please visit my ResearchGate for more details about ADP: https://www.researchgate.net/profile/Tao_Bian2

Features
-----------------

#### value iteration (VI) and policy iteration (PI) for continuous-time linear time invariant systems:

```c++
ALgorithmADP.h
```

#### VI and PI based ADP algorithms for continuous-time linear time invariant systems

```c++
ControllerADP.h
  ControllerVI.h
  ControllerPI.h
```

#### Different matrix classes and operations

```c++
Matrix.h
  SquareMatrix.h
    SymmetricMatrix.h
      Diagonal.h
MatrixCalc.h
```

#### Linear dynamical systems

```c++
Dynamical.h
```

How to use:
-----------------
A simple demo is provided in
```c++
main.cpp
```
An application of ADPlib in computer games can be found at: https://github.com/lyokka/ADP_lib/tree/master/ADPgame

#### Instruction for \*nix systems:
```bash
make all
./ADPsolver
```

Future work:
-----------------
1. discrete-time ADP
2. multi-thread implementation
3. nonlinear ADP
