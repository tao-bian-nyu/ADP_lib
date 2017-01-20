#ifndef DIAGNONAL_H
#define DIAGNONAL_H
#include <math.h>
#include "Diagonal.h"

namespace ADP
{

	Diagonal::Diagonal(const int ncol, const double val)
		:SymmetricMatrix(ncol)
	{
		for (unsigned int j=1; j<=ncol;j++)
			this->operator()(j,j) = val;

	}


	Diagonal::Diagonal(const Matrix& mat)
	:SymmetricMatrix((SquareMatrix(mat)).diag()){}
	
	Diagonal::Diagonal(const std::vector<double> & input)
		:SymmetricMatrix(input.size())
	{
		std::vector<double>::const_iterator iter = input.begin();

		for (unsigned int j=1; j<=input.size();j++)
		{
			this->operator()(j,j)=*iter;
			iter++;
		}

	}

}
#endif
