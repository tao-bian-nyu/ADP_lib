#ifndef ALGORITHMPI_H
#define ALGORITHMPI_H

#include <vector>
//#include "Matrix.h"
//#include "SymmetricMatrix.h"
//#include "SquareMatrix.h"
#include "Diagonal.h"
#include "AlgorithmADP.h"

namespace ADP
{

	class AlgorithmPI: public AlgorithmADP
	{
		public:

			AlgorithmPI(const SymmetricMatrix& Q, const SymmetricMatrix& R, const Matrix& K0);
			AlgorithmPI(): mQ(), mR(), mP(), mK0(), mK(), mk(1){};
			virtual const std::vector<Matrix>& offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N=20000, const double eps=1e-5);
			virtual const std::vector<Matrix>& online(const std::vector<double>& vec0, const std::vector<double>& vec1, const std::vector<double>& vec2, std::shared_ptr<Matrix> mBigr, SymmetricMatrix& mThetaInv, std::vector<double>& mBigV);
			virtual std::shared_ptr<AlgorithmADP> Creat(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, const Matrix& K0, Step* stepf, const double bound=10) const;
			virtual ~AlgorithmPI(){};
			virtual void disp() const;
		private:
			SymmetricMatrix mQ;
			SymmetricMatrix mR;
			SymmetricMatrix mP;
			Matrix mK0;
			Matrix mK;
			unsigned int mk;
			std::vector<Matrix> mResult;
			virtual void onlineI(const std::vector<double>&  vec);

	};
}
#endif
