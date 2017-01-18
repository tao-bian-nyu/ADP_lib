#ifndef ALGORITHMADP_H
#define ALGORITHMADP_H


#include <vector>
#include "Matrix.h"
#include "SymmetricMatrix.h"
#include "SquareMatrix.h"
#include "Diagonal.h"
#include "MatrixCalc.h"
#include "Step.h"

namespace ADP
{

	class AlgorithmADP
	{
		public:
			//virtual std::vector<Matrix> offline(const SquareMatrix& sysA, const Matrix& sysB, const SymmetricMatrix& Q, const SymmetricMatrix& R, double eps=0.001)=0;
			virtual std::vector<Matrix> offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N=20000, const double eps=0.001)=0;
			virtual std::vector<Matrix> onlineB(const std::vector<double>&  state, const std::vector<double>& input) = 0;
			virtual std::vector<Matrix> onlineI(const std::vector<double>&  vec)=0;
			virtual std::shared_ptr<AlgorithmADP> Creat(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, const Matrix& K0, Step* stepf, const double bound=10)=0;
			virtual void resetStep(){};
			//virtual void  online(const std::vector<double>&  state, const std::vector<double>& input) = 0;
			virtual ~AlgorithmADP(){} ;

	};


}

#endif
