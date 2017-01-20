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
			AlgorithmVI(): mQ(), mR(), mP0(), mP(), mbound(0), mk(1), mStep(nullptr){};
			virtual const std::vector<Matrix>& offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N=20000, const double eps=0.001);
			virtual const std::vector<Matrix>& online(const std::vector<double>& vec0, const std::vector<double>& vec1, const std::vector<double>& vec2, std::shared_ptr<Matrix> mBigr, SymmetricMatrix& mThetaInv, std::vector<double>& mBigV);
			virtual std::shared_ptr<AlgorithmADP> Creat(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, const Matrix& K0, Step* stepf, const double bound=10) const;
			virtual void resetStep();
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
			virtual void onlineI(const std::vector<double>&  vec);

	};
}
#endif
