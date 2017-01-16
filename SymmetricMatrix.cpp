//#include <vector>
#include <math.h>
#include <iostream>
#include "SymmetricMatrix.h"
//#include "SquareMatrix.h"
//#include "Matrix.h"
#include "MatrixCalc.h"

namespace ADP
{

	//SymmetricMatrix::SymmetricMatrix(const int ncol, const double val)
	//:SquareMatrix(ncol, val)
	//{
	//}

	//SymmetricMatrix::SymmetricMatrix(const SymmetricMatrix& mat)
	//:SquareMatrix(mat)
	//{
	//}

	SymmetricMatrix::SymmetricMatrix(const std::vector<double> & input)
		:SquareMatrix(symget_dim(input))
	{
		std::vector<double>::const_iterator iter = input.begin();
		const int ncol = symget_dim(input);

		for (int j=1; j<=ncol;j++)
		{
			this->operator()(j,j)=*iter;
			++iter;
			for (int i=j+1; i<=ncol; i++)
			{
				this->operator()(j,i)=*iter;
				this->operator()(i,j)=*iter;
				++iter;
			}
		}
	}

	SymmetricMatrix::SymmetricMatrix(const Matrix& mat)
		:SquareMatrix((SquareMatrix(mat)+SquareMatrix(mat.t()))*0.5)
		 //:SquareMatrix(mat)
	{
		//if(mat.size()[0]!=mat.size()[1])
		//std::cout <<"warning! not symmetric matrix" << std::endl;
	}
	//*this = *this + t(mat);
	//const std::vector<double> Matrix::vecs() const
	//{
	//std::vector<double> vecOut;
	//for(int i=1;i<=mnrow;i++)
	//{
	//vecOut.push_back(this->operator()(i,i));
	//for(int j=i+1;j<=mncol;j++)
	//{
	//vecOut.push_back(2*(this->operator()(i,j)));
	//}
	//}
	//return vecOut;
	//}

	const bool SymmetricMatrix::operator>(const double input)
	{
		try{
		Matrix mat=cholesky(*this-input);
		}
		catch(...){return 0;}
		return 1;
	}
	const std::vector<double> SymmetricMatrix::vecs() const
	{
		std::vector<double> vecout;

		for (int j=1; j<=this->size()[0];j++)
		{
			vecout.push_back(this->operator()(j,j));
			for (int i=j+1; i<=this->size()[0]; i++)
			{
				vecout.push_back(2*(this->operator()(j,i)));
			}
		}
		return vecout;
	}

	const int SymmetricMatrix::symget_dim(const std::vector<double> & input) const
	{
		const int n = floor((sqrt( 1 + 8*static_cast<double>(input.size()) ) - 1 ) / 2);
		if (n*(n+1)/2!=input.size())
			std::cout<<"warning! not match for symmetric matrix." << std::endl;
		return n; 
	}

}
