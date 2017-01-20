#ifndef ALGORITHMRLS_H
#define ALGORITHMRLS_H

#include "SquareMatrix.h"

namespace ADP
{

	class AlgorithmRLS
	{
		public:
			AlgorithmRLS(const Matrix& Phi, const std::vector<double>& d, const long double epsilon = 1e-10);
			const std::vector<double>& disp() const;
			void  print() const;
		private: 
			AlgorithmRLS(const AlgorithmRLS& );
			const AlgorithmRLS& operator=(const AlgorithmRLS& ); 
			std::vector<double> mW;
			SquareMatrix mP;
			long double mEpsilon;

	};


}

#endif
