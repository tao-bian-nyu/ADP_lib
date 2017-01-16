#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
//#include "SquareMatrix.h" 

namespace ADP
{
	//class SquareMatrix
	class Matrix
	{
		public:
			explicit Matrix(const unsigned int ncol=1,const unsigned int nrow=1, const double val=0.0);
			Matrix(const std::vector<double> & input, const unsigned int ncol=1); //column by column
			Matrix(const Matrix& mat):mncol(mat.mncol),mnrow(mat.mnrow),matrix(mat.matrix){};
			virtual const std::vector<int> size() const; // [0]=mncol, [1]=mnrow
			double& operator()(const unsigned int colIdx, const unsigned int rowIdx);
			const double& operator()(const unsigned int colIdx, const unsigned int rowIdx) const;
			virtual const Matrix& operator=(const Matrix& mat);
			void col(const unsigned int col_idx, const std::vector<double>& vec);
			void row(const unsigned int row_idx, const std::vector<double>& vec);
			const Matrix  operator-() const;
			const std::vector<double> col(const unsigned int col_idx) const;
			const std::vector<double> row(const unsigned int row_idx) const;
			const std::vector<double>& vec() const;
			//const std::vector<double> vecs() const;
			const Matrix t() const;
			const Matrix& add(const double val);
			const double F() const;
			void disp() const;
			//const Matrix inv() const;
			explicit operator double() const;
			//operator SquareMatrix() const;
			virtual ~Matrix(){};
		private:
			unsigned int mncol;
			unsigned int mnrow;
			std::vector<double> matrix;
	};

}
#endif
