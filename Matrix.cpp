//#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include "Matrix.h"
//#include "AlgorithmRLS.h"

namespace ADP
{


	Matrix::Matrix(const unsigned int ncol, const unsigned int nrow, const double val)
		:mncol(ncol),
		mnrow(nrow),
		matrix(ncol*nrow,0)
	{
		for (auto it=matrix.begin();it!=matrix.end();++it){*it=val;}
		
	}

	Matrix::Matrix(const std::vector<double> & input, const unsigned int ncol)
		:mncol(0),
		mnrow(0),
		matrix(input)
	{
		mncol = std::min(int(ncol),static_cast<int>(input.size()));
		mnrow = input.size()/mncol;

		if(mncol*mnrow!=input.size())
		{
			std::cout << "omit tail" << std::endl;
			matrix.resize(mncol*mnrow);
		}

	}


	const double& Matrix::operator()(const unsigned int colIdx, const unsigned int rowIdx) const
	{
		//return matrix[(ncol-1)*mnrow+nrow-1];
		if (colIdx>mnrow || colIdx<1)
		{
			std::cout <<"error column index" << std::endl;
			//throw invalid_argument("error column index");
			throw;
			//throw "error column index";
		}
		if (rowIdx>mncol || rowIdx<1)
		{
			std::cout <<"error row index" << std::endl;
			throw;
		}
		return matrix[(rowIdx-1)*mnrow+colIdx-1];
	}

	double& Matrix::operator()(const unsigned int colIdx, const unsigned int rowIdx)  // A_ncol,nrow
	{
		return	const_cast<double&>(
					static_cast<const Matrix&>(*this)(colIdx,rowIdx));
		//return const_cast<double&> (this->operator()(colIdx,rowIdx));
	}

	const Matrix& Matrix::operator=(const Matrix& mat)
	{
		if(this == &mat) return *this;

		mncol = mat.mncol;
		mnrow = mat.mnrow;
		matrix = mat.matrix;
		//const int kend = mncol*mnrow;
		//// matrix = new double[kend];
		//for (int k=0;k<kend;k++){matrix[k]=mat.matrix[k];}
		return *this;
	}



	void Matrix::col(const unsigned int col_idx, const std::vector<double>& vec)
	{
		if(vec.size()!=mnrow){
			std::cout <<"warning! unmatched length" << std::endl;
		}
		std::vector<double>::const_iterator it = vec.begin();
		for (unsigned int k=0;k<mnrow;k++)
		{
			matrix[col_idx*mnrow-mnrow+k]=*it;
			++it;
			if (it==vec.end()) break;
		}

	}

	const std::vector<double> Matrix::col(const unsigned int col_idx) const
	{
		std::vector<double> vec(mnrow,0.0);
		std::vector<double>::iterator it = vec.begin();
		for (unsigned int k=0;k<mnrow;k++)
		{
			*it = matrix[col_idx*mnrow-mnrow+k];
			++it;
		}
		return vec;

	}

	void Matrix::row(const unsigned int row_idx, const std::vector<double>& vec) 
	{
		//std::vector<double>::const_iterator it = vec.begin();
		if(vec.size()!=mncol){
			std::cout <<"warning! unmatched length" << std::endl;
		}

		std::vector<double>::const_iterator it = vec.begin();
		for (unsigned int k=0;k<mncol;k++)
		{
			matrix[k*mnrow+row_idx-1]=*it;
			++it;
			if (it==vec.end()) break;

		}

	}


	const std::vector<double> Matrix::row(const unsigned int row_idx) const
	{
		std::vector<double> vec(mncol,0.0);
		std::vector<double>::iterator it = vec.begin();
		for (unsigned int k=0;k<mncol;k++)
		{
			*it = matrix[k*mnrow+row_idx-1];
			++it;
		}
		return vec;

	}

	const Matrix Matrix::operator-() const
	{

		std::vector<double> vecOut(matrix);
		transform(vecOut.begin(), vecOut.end(), vecOut.begin(), [=](double i){return -i;});
		Matrix matOut(vecOut,mncol);
		return matOut;
	}

	const std::vector<double>& Matrix::vec() const
	{
		return matrix;
	}


	const Matrix Matrix::t() const
	{
		//int ncol = mnrow;
		//int nrow = mncol;

		Matrix matout(mnrow, mncol, 0.0);

		for (unsigned int i=1;i<=mncol;i++)
			matout.row(i, (this->col(i)));
		return matout;
	}



	//const double Matrix::F() const
	//{
		//double out = accumulate(matrix.begin(), matrix.end(), 0.0, [](double i, double j){return i + j*j;});
		//return sqrt(out);

	//}

	const Matrix& Matrix::add(const double val) 
	{
		transform(matrix.begin(), matrix.end(), matrix.begin(), [=](double i){return i+val;});
		return *this;
	}



	const std::vector<int>  Matrix::size() const
	{
		std::vector<int> out;
		out.push_back(mncol);
		out.push_back(mnrow);
		return out;
	}


	Matrix::operator double() const
	{
		if (mncol==1 && mnrow==1)
			return this->operator()(1,1);
		else{
			std::cout << "cannot transfer" << std::endl;
			throw;
		}
	}

	void Matrix::disp() const
	{
		std::cout << std::endl;
		for(unsigned int i=1;i<=mnrow;i++)
		{
			for(unsigned int j=1;j<=mncol;j++)
				std::cout << std::setw(5) << this->operator()(i,j) << ',';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}
