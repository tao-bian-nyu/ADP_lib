#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H
#include "Matrix.h"
namespace ADP
{

	class SquareMatrix: public Matrix
	{
		public:
			explicit SquareMatrix(const int ncol=1, const double val=0.0):Matrix(ncol, ncol, val){};
			SquareMatrix(const SquareMatrix& mat): Matrix(mat){};
			SquareMatrix(const Matrix& mat): Matrix(check_mat(mat)){};
			explicit SquareMatrix(const std::vector<double> & input);
			const SquareMatrix inv() const;                                                              // matrix inverse
			const std::vector<double> diag() const;                                                      // takeout diagonal as a vector
			virtual ~SquareMatrix(){};
		private:
			const int squareGet_dim(const std::vector<double> & input) const;
			const Matrix check_mat(const Matrix& mat) const;
	};
}
#endif
