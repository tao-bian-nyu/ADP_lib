#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
//#include "SquareMatrix.h" 

namespace ADP
{
	class Matrix
	{
		public:
			explicit Matrix(const unsigned int ncol=1, const unsigned int nrow=1, const double val=0.0);             // (# columns, # of rows, each element)
			Matrix(const std::vector<double> & input, const unsigned int ncol=1);                                    // convert vector to matrix column by column
			Matrix(const Matrix& mat):mncol(mat.mncol),mnrow(mat.mnrow),matrix(mat.matrix){};          
			const std::vector<unsigned int> size() const;                                                                     // [0]=mncol, [1]=mnrow
			double& operator()(const unsigned int colIdx, const unsigned int rowIdx);
			const double& operator()(const unsigned int colIdx, const unsigned int rowIdx) const;
			const Matrix& operator=(const Matrix& mat);
			void col(const unsigned int col_idx, const std::vector<double>& vec);                                    // put vec into col_idx-th column
			void row(const unsigned int row_idx, const std::vector<double>& vec);                                    // put vec into row_idx-th row 
			const Matrix  operator-() const;
			const std::vector<double> col(const unsigned int col_idx) const;
			const std::vector<double> row(const unsigned int row_idx) const;
			const std::vector<double>& vec() const;                                                                  // convert matrix into vector column by column
			const Matrix t() const;
			void clean();
			const Matrix& add(const double val);
			const Matrix& add(const Matrix mat);
			void disp() const;
			explicit operator double() const;
			virtual ~Matrix(){};
		private:
			unsigned int mncol;
			unsigned int mnrow;
			//std::vector<unsigned int> msize;
			//static std::vector<double> mcol;
			//static std::vector<double> mrow;
			std::vector<double> matrix;
	};

}
#endif
