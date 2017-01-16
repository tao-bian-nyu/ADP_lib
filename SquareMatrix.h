#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H
#include "Matrix.h"
namespace ADP
{

	class SquareMatrix: public Matrix
	{
		public:
			explicit SquareMatrix(const int ncol, const double val=0.0):Matrix(ncol, ncol, val){};
			SquareMatrix(const SquareMatrix& mat): Matrix(mat){};
			SquareMatrix(const Matrix& mat): Matrix(check_mat(mat)){};
			explicit SquareMatrix(const std::vector<double> & input);
			//virtual void add(const SquareMatrix& mat);
			//virtual void add(const double val);
			//virtual SquareMatrix& operator=(SquareMatrix& mat);
			virtual const SquareMatrix inv() const;
			virtual const std::vector<double> diag() const;
			virtual ~SquareMatrix(){};
		private:
			const int squareGet_dim(const std::vector<double> & input) const;
			const Matrix check_mat(const Matrix& mat) const;
	};
}
#endif
