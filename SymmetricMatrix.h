#ifndef SYMMETRICMATRIX_H
#define SYMMETRICMATRIX_H
//#include <vector>
#include "SquareMatrix.h"

namespace ADP
{

  class SymmetricMatrix: public SquareMatrix
  {
  public:
    SymmetricMatrix(const int ncol=1, const double val=0.0): SquareMatrix(ncol, val){};
    SymmetricMatrix(const SymmetricMatrix& mat):SquareMatrix(mat){};
    SymmetricMatrix(const Matrix& mat);
    SymmetricMatrix(const std::vector<double> & input);
    const bool operator>(const double input);                                                     // check symmetric positive definite matrix
    const std::vector<double> vecs() const;                                                       // convert lower triangle into a vector
    virtual ~SymmetricMatrix(){};

  private:
    const int symget_dim(const std::vector<double> & input) const;
  };

}
#endif
