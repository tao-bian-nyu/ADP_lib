#include "MatrixCalc.h"
#include "Matrix.h"
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>


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
		SquareMatrix matout(lhs_mat);
		// matout.matrix = new double[kend];
		const unsigned int n = matout.size()[0];
		Diagonal rhs_mat(n,rhs_num);
		matout = matout + rhs_mat;

		//for (int i=1;i<=ncol;i++){
			//// matout.matrix[k]+=rhsmat.matrix[k];
			//matout(i,i)+=rhs_mat(i,i);
		//}
		return matout;
	}

	const Matrix operator+(const Matrix& lhs_mat, const Matrix& rhs_mat)
	{
		if (lhs_mat.size()[0]==rhs_mat.size()[0] && lhs_mat.size()[1]==rhs_mat.size()[1])
		{
			Matrix matout(lhs_mat);
			for (unsigned int i=1;i<=matout.size()[1];i++)
				for (unsigned int j=1;j<=matout.size()[0];j++)
					matout(i,j)+=rhs_mat(i,j);
			return matout;

		}

		std::cout << "error! two matrices are not match" << std::endl;
		throw;
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
			//return lhs_vec;
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
			//return lhs_vec;
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
			//return 0;
		}

		return std::inner_product(lhs_vec.begin(), lhs_vec.end(), rhs_vec.begin(),0.0);

		//std::vector<double> vecOut(lhs_vec);
		//std::transform(vecOut.begin(), vecOut.end(), rhs_vec.begin(), vecOut.begin(), std::multiplies<double>());
		////double out = std::accumulate(vecOut.begin(),vecOut.end(),0.0);
		//return std::accumulate(vecOut.begin(),vecOut.end(),0.0);

	}

	const Matrix operator*(const Matrix& lhs_mat, const double rhs_num)
	{

		Matrix matout(lhs_mat);

		for (unsigned int i=1;i<=matout.size()[1];i++)
			for (unsigned int j=1;j<=matout.size()[0];j++)
				matout(i,j)*=rhs_num;

		return matout;
	}

	const Matrix operator*(const Matrix& lhs_mat, const Matrix& rhs_mat)
	{
		//if (rhs_mat.size()[1]!=lhs_mat.size()[0]) return rhs_mat;
		if (rhs_mat.size()[1]!=lhs_mat.size()[0])
		{
			std::cout <<"error! matrices do not match" << std::endl;
			throw;
		}

		Matrix matout(rhs_mat.size()[0],lhs_mat.size()[1]);
		for (unsigned int i=1;i<=matout.size()[1];i++)
			for (unsigned int j=1;j<=matout.size()[0];j++)
				matout(i,j) = lhs_mat.row(i)*rhs_mat.col(j);
		return matout;
	}

	const Matrix operator*(const double lhs_num, const Matrix& rhs_mat)
	{
		return rhs_mat * lhs_num;

	}

	const Matrix prod(const std::vector<double> & input1, const std::vector<double> & input2)
	{
		Matrix matOut(input2.size(),input1.size());
		for (unsigned int i=1;i<=input1.size();i++)
			for (unsigned int j=1;j<=input2.size();j++)
				matOut(i,j) = input1[i-1]*input2[j-1];

		return matOut;

	}

	const Matrix kProd(const Matrix&  mat1, const Matrix&  mat2)
	{
		Matrix matOut(mat1.size()[0]*mat2.size()[0], mat1.size()[1]*mat2.size()[1]);

		for (unsigned int i=1;i<=mat1.size()[1];i++)
			for (unsigned int j=1;j<=mat1.size()[0];j++)
				for (unsigned int k=1;k<=mat2.size()[1];k++)
					for (unsigned int n=1;n<=mat2.size()[0];n++)
					{
						const unsigned int left = (i-1)*mat2.size()[1]+k;
						const unsigned int right = (j-1)*mat2.size()[0]+n;
						matOut(left,right) = mat1(i,j)*mat2(k,n);
					}

		return matOut;

	}


	const SquareMatrix kSum(const SquareMatrix&  mat1, const SquareMatrix&  mat2)
	{
		Diagonal I1(mat2.size()[0]);
		Diagonal I2(mat1.size()[0]);

		return kProd(mat1,I1)+kProd(I2,mat2);
	}



	const Matrix t(const Matrix& mat)
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
		auto it = vec.begin();

		std::cout<<std::endl;
		for (;it!=vec.end();++it)
		{
			std::cout<< std::setw(5) << *it << ',';
		}
		std::cout<<std::endl;

	}
	const Matrix inv(const SquareMatrix& mat)
	{
		return mat.inv();
	}

	const Matrix cholesky(const SymmetricMatrix&  mat)
	{
		unsigned int n = mat.size()[0];
		Matrix matOut(n,n,0);

		for (unsigned int i=1;i<=n;++i)
			for (unsigned int j=1;j<=i;++j){
				double s = matOut.row(i) * matOut.row(j); 
				matOut(i,j) = (i==j)?
					sqrt(mat(i,i)-s):
					1.0/matOut(j,j)*(mat(i,j)-s);
				if (std::isfinite(matOut(i,j))==0)
				{
					//std::cout << "not semipositive definite!" << std::endl;
					throw std::invalid_argument("not semipositive definite!");
				}
				//std::cout<<std::isfinite(matOut(i,j));
			}
		return matOut;
	}
	//void matDisplay(const Matrix&  mat)
	//{
		//for(int i=1;i<=mat.size()[1];i++)
		//{
			//for(int j=1;j<=mat.size()[0];j++)
				//std::cout << std::setw(5) << mat(i,j) << ',';
			//std::cout << std::endl;
		//}
	//}
}
