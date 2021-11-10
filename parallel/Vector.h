/*!
 * @file Vector.h
 *
 * @brief Defines a Vector class.
 *
 * @authors Gabriel Suau, Remi Pegouret, Geoffrey Lebaud
 *
 * @version 0.1.0
 *
 * @copyright © 2021 Gabriel Suau
 * @copyright © 2021 Remi Pegouret
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


#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>

/*!
 * @class DVector
 *
 * @brief Represents a vector of doubles.
 *
 * @details This class is derived from the std::vector<double> class of the STL.
 * We added method and operators to perform basic operations such as : addition, 
 * substraction, multiplication with a scalar, and dot product.
 */
class DVector: public std::vector<double>
{
public:
  /*! @brief Empty constructor. */
  DVector();
  
  /*! @brief Constructs a vector full of zeros of size count. */
  DVector(size_type count);

  // Additional methods
  /*! @brief Adds the DVector vec to this DVector. */
  DVector add(const DVector& vec);
  
  /*! @brief Subtracts the DVector vec to this DVector. */
  DVector sub(const DVector& vec);
  
  /*! @brief Computes the dot product between the DVector vec and this DVector. */
  double dot(const DVector& vec);

  // Print the vector (for debugging purpose)
  /*! @brief Prints this DVector (only used for debugging purposes). */
  void print() const;
};

/*!
 * @brief Overloads the << operator for DVectors.
 */
std::ostream& operator<< (std::ostream &os, const DVector& v);
/*!
 * @brief Computes the sum of two DVectors.
 */
DVector operator+ (const DVector& u, const DVector& v);
/*!
 * @brief Computes the difference between two DVectors.
 */
DVector operator- (const DVector& u, const DVector& v);
/*!
 * @brief Computes the scalar-vector product.
 */
DVector operator* (double alpha, const DVector& u);
/*!
 * @brief Computes the vector-scalar product.
 */
DVector operator* (const DVector& u, double alpha);

#endif //VECTOR_H
