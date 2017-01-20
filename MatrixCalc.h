#ifndef MATRIXCALC_H 
#define MATRIXCALC_H 


//#include "MatrixCalc.h"
//#include "Matrix.h"
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
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
	const Matrix operator/(const Matrix& lhs_mat, const double rhs_num);
	const Matrix operator*(const double lhs_num, const Matrix& rhs_mat);

	const Matrix prod(const std::vector<double> & input1, const std::vector<double> & input2);
	const Matrix kProd(const Matrix&  mat1, const Matrix&  mat2);
	const SquareMatrix kSum(const SquareMatrix&  mat1, const SquareMatrix&  mat2);
	const Matrix T(const Matrix&  mat);
	const std::vector<double> vec(const Matrix&  mat);
	void disp(const std::vector<double>& vec);
	void disp(const Matrix& mat);

	const Matrix inv(const SquareMatrix& mat);
	const std::vector<double> diag(const SquareMatrix& mat);
	const std::vector<double> vecs(const SymmetricMatrix&  mat);
	const SquareMatrix cholesky(const SymmetricMatrix&  mat);

	template<char E = 'E'>
	const double norm(const Matrix& mat);                                                                     // Euclidean norm
	
	template<>
	const double norm<'F'>(const Matrix& mat);                                                                // Frobenius norm

	template<>
	const double norm<'I'>(const Matrix& mat);                                                                // Infinite norm

	template<>
	const double norm<'1'>(const Matrix& mat);                                                                // 1norm

	template<int P = 2>
	const double norm(const std::vector<double>& vec)                                                         // p norm
	{
		double sum = 0;
		std::for_each(vec.begin(), vec.end(), [&](double i) { sum += pow(std::abs(i),P);} );
		return pow(sum,1/P);
	}


}

//#include "MatrixCalc.cpp"
#endif
