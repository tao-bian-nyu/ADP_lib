#ifndef ALGORITHMADP_H
#define ALGORITHMADP_H


#include <vector>
#include <memory>
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
			virtual const std::vector<Matrix>& offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N=20000, const double eps=0.001)=0;
			virtual const std::vector<Matrix>& online(const std::vector<double>& vec0, const std::vector<double>& vec1, const std::vector<double>& vec2, std::shared_ptr<Matrix> mBigr, SymmetricMatrix& mThetaInv, std::vector<double>& mBigV)=0;
			virtual std::shared_ptr<AlgorithmADP> Creat(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, const Matrix& K0, Step* stepf, const double bound=10) const =0;
			virtual void resetStep(){};
			virtual ~AlgorithmADP(){};
			virtual void disp() const=0;
		private:
			virtual void onlineI(const std::vector<double>&  vec)=0;

	};


}

#endif
