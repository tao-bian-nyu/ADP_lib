//#include <vector>
//#include "Matrix.h"
//#include "SymmetricMatrix.h"
//#include "SquareMatrix.h"
#ifndef MATRIXCALC_H 
#define MATRIXCALC_H 


#include "Diagonal.h"

namespace ADP
{

	const std::vector<double> operator+(const std::vector<double>& lhs_vec, const double rhs_num);
	const std::vector<double> operator-(const std::vector<double>& lhs_vec, const double rhs_num);
	const std::vector<double> operator*(const std::vector<double>& lhs_vec, const double rhs_num);
	const std::vector<double> operator/(const std::vector<double>& lhs_vec, const double rhs_num);
	const std::vector<double> operator+(const double lhs_num, const std::vector<double>& rhs_vec);
	const std::vector<double> operator-(const double lhs_num, const std::vector<double>& rhs_vec);
	const std::vector<double> operator*(const double lhs_num, const std::vector<double>& rhs_vec);
	const std::vector<double> operator+(const std::vector<double>& lhs_vec, const std::vector<double>& rhs_vec);
	const std::vector<double> operator-(const std::vector<double>& lhs_vec, const std::vector<double>& rhs_vec);
	const double operator*(const std::vector<double>& lhs_vec, const std::vector<double>& rhs_vec);

	const Matrix operator+(const Matrix& lhs_mat, const Matrix& rhs_mat);
	const SquareMatrix operator+(const double lhs_num, const SquareMatrix& rhs_mat);
	const SquareMatrix operator+(const SquareMatrix& lhs_mat, const double rhs_num);
	const Matrix add(const double lhs_num, const Matrix& rhs_mat);
	const Matrix add(const Matrix& lhs_mat, const double rhs_num);
	const Matrix add(const Matrix& lhs_mat, const Matrix& rhs_mat);

	const Matrix operator-(const Matrix& lhs_mat, const Matrix& rhs_mat);
	const SquareMatrix operator-(const double lhs_num, const SquareMatrix& rhs_mat);
	const SquareMatrix operator-(const SquareMatrix& lhs_mat, const double rhs_num);

	const Matrix operator*(const Matrix& lhs_mat, const Matrix& rhs_mat);
	const Matrix operator*(const Matrix& lhs_mat, const double rhs_num);
	const Matrix operator*(const double lhs_num, const Matrix& rhs_mat);

	const Matrix prod(const std::vector<double> & input1, const std::vector<double> & input2);
	const Matrix kProd(const Matrix&  mat1, const Matrix&  mat2);
	const SquareMatrix kSum(const SquareMatrix&  mat1, const SquareMatrix&  mat2);
	const Matrix t(const Matrix&  mat);
	const std::vector<double> vec(const Matrix&  mat);
	void disp(const std::vector<double>& vec);
	void disp(const Matrix& mat);

	const Matrix inv(const SquareMatrix& mat);
	const std::vector<double> diag(const SquareMatrix& mat);
	const std::vector<double> vecs(const SymmetricMatrix&  mat);
	const Matrix cholesky(const SymmetricMatrix&  mat);
}

#endif
