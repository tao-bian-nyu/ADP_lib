#include <vector>
#include <math.h>
//#include <iostream>
#include "Controllers.h"
#include "Matrix.h"
#include "MatrixCalc.h"
#include "AlgorithmADP.h"
#include "Step.h"
#include "Others.h"
//#include "AlgorithmVI.h"
//#include "AlgorithmPI.h"

namespace ADP{

	const std::vector<double> sinusoidal(const unsigned int m, const double t)
	{
		return std::vector<double>(m,sin(t));
	}

	const std::vector<double> linInput(const std::vector<double>& x, const int m, const double t){
		Matrix K(x.size(),m,1);
		std::vector<double> u = (K * x).vec();
		return u;
	}
}
