/*!
 * @file Laplacian.h
 *
 * @brief Defines a class representing the discrete laplacian matrix.
 *
 * @authors Gabriel Suau, Lucas Trautmann, Geoffrey Lebaud
 *
 * @version 0.1.0
 *
 * @copyright © 2021 Gabriel Suau
 * @copyright © 2021 Lucas Trautmann
 * @copyright © 2021 Geoffrey Lebaud
 *
 * @copyright This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "MPIUtils.h"
#include "DataFile.h"
#include "Function.h"
#include "Vector.h"
#include <fstream>

/*!
 * @class Laplacian
 *
 * @brief Represents the discrete 2D Laplacian matrix
 *
 * @details A Laplacian object contains everything that is needed to create the discrete Laplacian matrix,
 * perform a matrix-vector product and solve a linear system involving this matrix with the
 * conjugate gradient method.
 *
 * @details This matrix is composed of tridiagonal blocks on its diagonal and diagonal blocks on its sub/super-diagonal.
 * This matrix is sparse, so it is stored in CSR format, and the algorithms for the matvec products and the conjugate gradient
 * are adapted.
 */
class Laplacian
{
private:
  DataFile* _DF; ///< Pointer to a DataFile object.
  Function* _function; ///< Pointer to a Function object.

  DVector AA, IA, JA; ///< CSR storage.
  int _NNZ; ///< Number of non zero elements in the matrix.
  int _Nx, _Ny; ///< Number of unknowns in the x and y directions.

public:
  /*! @brief Constructs an empty Laplacian object. */
  Laplacian();

  /*! @brief Constructs a Laplacian object using a DataFile object and a Function object. */
  Laplacian(DataFile* DF, Function* function);

  /*! @brief Default destructor. */
  ~Laplacian() = default;

  /*! @brief Initializes an already constructed Laplacian object. */
  void Initialize();

  /*! @brief Initializes an empty constructed Laplacian object. */
  void Initialize(DataFile* DF, Function* function);

  // Getters
  int getNx() const {return _Nx;}
  int getNy() const {return _Ny;};

  /*! @brief Performs a matrix vector product. */
  DVector matVecProd(const DVector& x);

  /*! @brief Solves the linear system Ax = b using the conjugate gradient method. */
  DVector solveConjGrad(const DVector& b, const DVector& x0, double tolerance, int maxIterations, std::ofstream& resFile);
};

#endif // LAPLACIAN_H
