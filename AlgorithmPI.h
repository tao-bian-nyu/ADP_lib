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

			AlgorithmPI(const SymmetricMatrix& Q, const SymmetricMatrix& R, const Matrix& K0): mQ(Q), mR(R), mP(Q), mK0(K0), mK(K0), mk(1){};
			AlgorithmPI(): mQ(), mR(), mP(), mK0(), mK(), mk(1){};
			virtual std::shared_ptr<AlgorithmADP> Creat(const SymmetricMatrix& Q, const SymmetricMatrix& R, const SymmetricMatrix& P0, const Matrix& K0, Step* stepf, const double bound=10);
			virtual std::vector<Matrix> offline(const SquareMatrix& sysA, const Matrix& sysB, const unsigned int N=20000, const double eps=1e-5);
			virtual std::vector<Matrix> onlineI(const std::vector<double>&  vec);
			virtual std::vector<Matrix> onlineB(const std::vector<double>&  state, const std::vector<double>& input){return mResult;};
			//virtual std::vector<Matrix> onlineI(const std::vector<double>&  state, const std::vector<double>& input){return mResult;};
			//virtual void online(const std::vector<double>&  state, const std::vector<double>& input){};
			virtual ~AlgorithmPI();
		private:
			SymmetricMatrix mQ;
			SymmetricMatrix mR;
			SymmetricMatrix mP;
			Matrix mK0;
			Matrix mK;
			unsigned int mk;
			std::vector<Matrix> mResult;

	};
}
#endif
