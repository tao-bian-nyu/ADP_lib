#ifndef RK_H
#define RK_H

#include "SquareMatrix.h"
#include "Controllers.h"

namespace ADP{
	class RK
	{
		public:
			typedef const std::vector<double> (*inputfun) (const std::vector<double>&, const int, const double);
			RK(const SquareMatrix A, const Matrix B, const double dt):mA(A),mB(B),mdt(dt){};
			const std::vector<double> linear(const std::vector<double>& x, Controllers* controller = nullptr, inputfun input = nullptr, const double t=0);
		private:
			SquareMatrix mA;
			Matrix mB;
			double mdt;
	};
}

#endif
