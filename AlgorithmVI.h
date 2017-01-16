#ifndef ALGORITHMVI_H
#define ALGORITHMVI_H

//#include <vector>
//#include "Matrix.h"
//#include "SymmetricMatrix.h"
//#include "SquareMatrix.h"
#include "Diagonal.h"
//#include "MatrixCalc.h"
#include "AlgorithmADP.h"
#include "Step.h"

namespace ADP
{

	class AlgorithmVI: public AlgorithmADP
	{
		public:

			AlgorithmVI(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, Step* stepf, const double bound=10);
			//AlgorithmVI(): mP0(0), mbound(0), mResult({0, 0, 0}){};
			virtual std::vector<Matrix> offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N=20000, const double eps=0.001);
			//virtual void offline(const SquareMatrix& sysA, const Matrix& sysB, const SquareMatrix& Q, const SquareMatrix& R, double eps=0.001);
			virtual std::vector<Matrix> onlineB(const std::vector<double>&  state, const std::vector<double>& input);
			virtual std::vector<Matrix> onlineI(const std::vector<double>&  vec);
			virtual void resetStep();
			//virtual void online(const std::vector<double>&  state, const std::vector<double>& input){};
			virtual ~AlgorithmVI(){};
		private:
			SymmetricMatrix mQ;
			SymmetricMatrix mR;
			SymmetricMatrix mP0;
			SymmetricMatrix mP;
			double mbound;
			unsigned int mk;
			std::vector<Matrix> mResult;
			Step* mStep;

	};
}
#endif
