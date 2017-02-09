#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include "MatrixCalc.h"
#include "Matrix.h"


namespace ADP
{

	const Matrix add(const Matrix& lhs_mat, const double rhs_num)
	{
		Matrix matOut(lhs_mat);
		matOut.add(rhs_num);
		return matOut;
	}


	const Matrix add(const double lhs_num, const Matrix& rhs_mat)
	{
		return add(rhs_mat, lhs_num);
	}

	const Matrix add(const Matrix& lhs_mat, const Matrix& rhs_mat)
	{
		return lhs_mat + rhs_mat;
	}


	
	const SquareMatrix operator+(const SquareMatrix& lhs_mat, const double rhs_num)
	{
		const unsigned int n = lhs_mat.size()[0];
		Diagonal rhs_mat(n,rhs_num);
		return lhs_mat + rhs_mat;
	}

	const Matrix operator+(const Matrix& lhs_mat, const Matrix& rhs_mat)
	{
		Matrix matout(lhs_mat);
		return matout.add(rhs_mat);
		//if (lhs_mat.size()[0]==rhs_mat.size()[0] && lhs_mat.size()[1]==rhs_mat.size()[1])
		//{
			//Matrix matout(lhs_mat);
			//for (unsigned int i=1;i<=matout.size()[1];i++)
				//for (unsigned int j=1;j<=matout.size()[0];j++)
					//matout(i,j)+=rhs_mat(i,j);
			//return matout;

		//}

		//std::cout << "error! two matrices are not match" << std::endl;
		//throw;
	}

	const SquareMatrix operator+(const double lhs_num, const SquareMatrix& rhs_mat)
	{
		return rhs_mat + lhs_num;

	}


	const SquareMatrix operator-(const SquareMatrix& lhs_mat, const double rhs_num)
	{
		return lhs_mat + (-rhs_num);
	}


	const SquareMatrix operator-(const double lhs_num, const SquareMatrix& rhs_mat)
	{
		return -(rhs_mat - lhs_num);

	}

	const Matrix operator-(const Matrix& lhs_mat, const Matrix& rhs_mat)
	{
		return lhs_mat + (-rhs_mat);
	}

	const std::vector<double> operator+(const std::vector<double>& lhs_vec, const double rhs_num)
	{
		std::vector<double> vecOut(lhs_vec);
		transform(vecOut.begin(), vecOut.end(), vecOut.begin(), [=](double i){return i+rhs_num;});
		return vecOut;
	}

	const std::vector<double> operator-(const std::vector<double>& lhs_vec, const double rhs_num)
	{
		std::vector<double> vecOut(lhs_vec);
		transform(vecOut.begin(), vecOut.end(), vecOut.begin(), [=](double i){return i-rhs_num;});
		return vecOut;
	}

	const std::vector<double> operator*(const std::vector<double>& lhs_vec, const double rhs_num)
	{
		std::vector<double> vecOut(lhs_vec);
		transform(vecOut.begin(), vecOut.end(), vecOut.begin(), [=](double i){return i*rhs_num;});
		return vecOut;
	}

	const std::vector<double> operator/(const std::vector<double>& lhs_vec, const double rhs_num)
	{
		std::vector<double> vecOut(lhs_vec);
		transform(vecOut.begin(), vecOut.end(), vecOut.begin(), [=](double i){return i/rhs_num;});
		return vecOut;
	}

	const std::vector<double> operator+(const double lhs_num, const std::vector<double>& rhs_vec)
	{
		return rhs_vec + lhs_num;
	}
	const std::vector<double> operator-(const double lhs_num, const std::vector<double>& rhs_vec)
	{
		return (rhs_vec * (-1)+ lhs_num);
	}
	const std::vector<double> operator*(const double lhs_num, const std::vector<double>& rhs_vec)
	{
		return rhs_vec * lhs_num;
	}

	const std::vector<double> operator+(const std::vector<double>& lhs_vec, const std::vector<double>& rhs_vec)
	{
		if (rhs_vec.size()!=lhs_vec.size()) 
		{
			std::cout <<"error! different vector length" << std::endl;
			throw;
		}
		std::vector<double> vecOut(lhs_vec);
		transform(vecOut.cbegin(), vecOut.cend(), rhs_vec.begin(), vecOut.begin(), std::plus<double>());
		return vecOut;
	}

	const std::vector<double> operator-(const std::vector<double>& lhs_vec, const std::vector<double>& rhs_vec)
	{
		if (rhs_vec.size()!=lhs_vec.size()) 
		{
			std::cout <<"error! different vector length" << std::endl;
			throw;
		}
		std::vector<double> vecOut(lhs_vec);
		std::transform(vecOut.begin(), vecOut.end(), rhs_vec.begin(), vecOut.begin(), std::minus<double>());

		return vecOut;

	}

	const double operator*(const std::vector<double>& lhs_vec, const std::vector<double>& rhs_vec)
	{
		if (rhs_vec.size()!=lhs_vec.size()) 
		{
			std::cout <<"error! different vector length" << std::endl;
			throw;
		}

		return std::inner_product(lhs_vec.begin(), lhs_vec.end(), rhs_vec.begin(),0.0);

	}

	const Matrix operator*(const Matrix& lhs_mat, const double rhs_num)
	{

		Matrix matout(lhs_mat);

		for (unsigned int i=1;i<=matout.size()[1];i++)
			for (unsigned int j=1;j<=matout.size()[0];j++)
				matout(i,j)*=rhs_num;

		return matout;
	}


	const Matrix operator/(const Matrix& lhs_mat, const double rhs_num)
	{
		return 1/rhs_num * lhs_mat;
	}



	const Matrix operator*(const Matrix& lhs_mat, const Matrix& rhs_mat)
	{
		auto rhsSize = rhs_mat.size();
		auto lhsSize = lhs_mat.size();
		if (rhsSize[1]!=lhsSize[0])
		{
			std::cout <<"error! matrices do not match" << std::endl;
			throw;
		}

		Matrix matout(rhsSize[0],lhsSize[1]);
		for (unsigned int i=1;i<=lhsSize[1];i++)
			for (unsigned int j=1;j<=rhsSize[0];j++)
				for (unsigned int k=1;k<=lhsSize[0];k++)
				//matout(i,j) = lhs_mat.row(i)*rhs_mat.col(j);
					matout(i,j) += lhs_mat(i,k) * rhs_mat(k,j);
		return matout;
	}

	const Matrix operator*(const double lhs_num, const Matrix& rhs_mat)
	{
		return rhs_mat * lhs_num;

	}

	const Matrix prod(const std::vector<double> & input1, const std::vector<double> & input2)
	{
		auto size1 = input1.size();
		auto size2 = input2.size();
		Matrix matOut(size2,size1);
		for (unsigned int i=1;i<=size1;i++)
			for (unsigned int j=1;j<=size2;j++)
				matOut(i,j) = input1[i-1]*input2[j-1];

		return matOut;

	}

	const Matrix kProd(const Matrix&  mat1, const Matrix&  mat2)
	{
		auto size1 = mat1.size();
		auto size2 = mat2.size();
		Matrix matOut(size1[0]*size2[0], size1[1]*size2[1]);

		for (unsigned int i=1;i<=size1[1];i++)
			for (unsigned int j=1;j<=size1[0];j++)
				for (unsigned int k=1;k<=size2[1];k++)
					for (unsigned int n=1;n<=size2[0];n++)
						matOut((i-1)*size2[1]+k,(j-1)*size2[0]+n) = mat1(i,j)*mat2(k,n);
		return matOut;

	}


	const SquareMatrix kSum(const SquareMatrix&  mat1, const SquareMatrix&  mat2)
	{
		Diagonal I1(mat2.size()[0]);
		Diagonal I2(mat1.size()[0]);

		return kProd(mat1,I1)+kProd(I2,mat2);
	}



	const Matrix T(const Matrix& mat)
	{
		return mat.t();
	}

	const std::vector<double> vec(const Matrix&  mat)
	{
		return mat.vec();
	}


	const std::vector<double> diag(const SquareMatrix& mat)
	{
		return mat.diag();
	}

	const std::vector<double> vecs(const SymmetricMatrix&  mat)
	{
		return mat.vecs();
	}

	void disp(const Matrix& mat)
	{
		mat.disp();
	}

	void disp(const std::vector<double>& vec)
	{
		std::cout<<std::endl;
		for (auto it = vec.begin();it!=vec.end();++it)
		{
			std::cout<< std::setw(5) << *it << ',';
		}
		std::cout<<std::endl;

	}
	const Matrix inv(const SquareMatrix& mat)
	{
		return mat.inv();
	}

	void clean(std::vector<double>& vec)
	{
		transform(vec.begin(), vec.end(), vec.begin(), [=](double i){return 0;});
	}

	void clean(Matrix& mat)
	{
		mat.clean();
	}

	const SquareMatrix cholesky(const SymmetricMatrix&  mat)
	{
		unsigned int n = mat.size()[0];
		SquareMatrix matOut(n,0);

		for (unsigned int i=1;i<=n;++i)
			for (unsigned int j=1;j<=i;++j){
				double s = matOut.row(i) * matOut.row(j); 
				matOut(i,j) = (i==j)?
					sqrt(mat(i,i)-s):
					1.0/matOut(j,j)*(mat(i,j)-s);
				if (std::isfinite(matOut(i,j))==0)
				{
					throw std::invalid_argument("not semipositive definite!");
				}
			}
		return matOut;
	}


	template<>
	const double norm<'E'>(const Matrix& mat) 
	{
		SquareMatrix choMat(cholesky(mat*T(mat)));
		std::vector<double> vec(diag(choMat));
		return *std::max_element(vec.begin(), vec.end());

	}


	template<>
	const double norm<'F'>(const Matrix& mat) 
	{
		std::vector<double> vecOut(vec(mat));
		double out = std::inner_product(vecOut.begin(), vecOut.end(), vecOut.begin(),0.0);
		return sqrt(out);

	}

	template<>
	const double norm<'I'>(const Matrix& mat) 
	{
		double out =0;
		std::vector<double> vecOut;

		for (unsigned int i=1;i<=mat.size()[1]; ++i)
		{
			vecOut = mat.row(i);
			double tem = std::inner_product(vecOut.begin(), vecOut.end(), vecOut.begin(),0.0);
			out = std::max(tem,out);
		}
		return out;

	}

	template<>
	const double norm<'1'>(const Matrix& mat) 
	{
		double out =0;
		std::vector<double> vecOut;

		for (unsigned int i=1;i<=mat.size()[0]; ++i)
		{
			vecOut = mat.col(i);
			double tem = std::inner_product(vecOut.begin(), vecOut.end(), vecOut.begin(),0.0);
			out = std::max(tem,out);
		}
		return out;

	}


}
