#ifndef OTHERS_H 
#define OTHERS_H 
#include <vector>
#include <math.h>
//#include <iostream>
#include "Controllers.h"
#include "Matrix.h"
#include "MatrixCalc.h"
#include "AlgorithmADP.h"
#include "Step.h"
//#include "AlgorithmVI.h"
//#include "AlgorithmPI.h"

namespace ADP{

	const std::vector<double> sinusoidal(const unsigned int m, const double t);
	const std::vector<double> linInput(const std::vector<double>& x, const int m, const double t=0);

}

#endif
