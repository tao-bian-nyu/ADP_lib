//#include <vector>
#include <math.h>
#include <iostream>
#include "SymmetricMatrix.h"
//#include "SquareMatrix.h"
//#include "Matrix.h"
#include "MatrixCalc.h"

namespace ADP
{


	SymmetricMatrix::SymmetricMatrix(const std::vector<double> & input)
		:SquareMatrix(symget_dim(input))
	{
		std::vector<double>::const_iterator iter = input.begin();
		const int ncol = symget_dim(input);

		for (unsigned int j=1; j<=ncol;j++)
		{
			this->operator()(j,j)=*iter;
			++iter;
			for (unsigned int i=j+1; i<=ncol; i++)
			{
				this->operator()(j,i)=*iter;
				this->operator()(i,j)=*iter;
				++iter;
			}
		}
	}

	SymmetricMatrix::SymmetricMatrix(const Matrix& mat)
		:SquareMatrix((SquareMatrix(mat)+SquareMatrix(T(mat)))*0.5){}

	const bool SymmetricMatrix::operator>(const double input)
	{
		try{
			cholesky(*this-input);
		}
		catch(...){return 0;}
		return 1;
	}
	const std::vector<double> SymmetricMatrix::vecs() const
	{
		std::vector<double> vecout;
		auto thisSize = this->size()[0];

		for (unsigned int j=1; j<=thisSize;j++)
		{
			vecout.push_back(this->operator()(j,j));
			for (unsigned int i=j+1; i<=thisSize; i++)
				vecout.push_back(2*(this->operator()(j,i)));
		}
		return vecout;
	}

	const int SymmetricMatrix::symget_dim(const std::vector<double> & input) const
	{
		const unsigned int n = floor((sqrt( 1 + 8*static_cast<double>(input.size()) ) - 1 ) / 2);

		if (n*(n+1)/2!=input.size())
			std::cout<<"warning! not match for symmetric matrix." << std::endl;
		return n; 
	}

}
