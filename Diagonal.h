#ifndef DIAGONAL_H
#define DIAGONAL_H
//#include <vector>
#include "SymmetricMatrix.h"

namespace ADP
{

  class Diagonal:public SymmetricMatrix
  {
  public:
    explicit Diagonal(const int ncol, const double val=1.0);
    Diagonal(const Diagonal& mat):SymmetricMatrix(mat){};
    Diagonal(const Matrix& mat);
    explicit Diagonal(const std::vector<double> & input);
    //virtual const int size() const;
    //virtual const SymmetricMatrix& operator=(const SymmetricMatrix& mat);

  //private:
    //const int symget_dim(const std::vector<double> & input);
  };

}
#endif
