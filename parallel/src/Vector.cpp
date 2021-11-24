/*!
 * @file Vector.cpp
 *
 * @brief Defines a Vector class.
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


#include "Vector.h"
#include <vector>
#include <iostream>

// Constructors
DVector::DVector():
  std::vector<double>()
{
}


// Constructors
DVector::DVector(size_type count):
  std::vector<double>(count)
{
}


// Additional methods
/*!
 * @param [in] vec Dvector that is added to this Dvector.
 * @par Modifies
 * this DVector
 * @return the new value of this DVector.
 */
DVector DVector::add(const DVector& vec)
{
  int i, N;

  N = vec.size();

  for (i = 0 ; i < N ; ++i) {
    this->operator[](i) += vec[i];
  }

  return *this;
}


/*!
 * @param [in] vec Dvector that is subtracted to this Dvector.
 * @par Modifies
 * this DVector
 * @return the new value of this DVector.
 */
DVector DVector::sub(const DVector& vec)
{
  int i, N;

  N = vec.size();

  for (i = 0 ; i < N ; ++i) {
    this->operator[](i) -= vec[i];
  }

  return *this;
}


/*!
 * @param [in] vec Dvector that is dotted to this Dvector.
 * @return the value of the dot product.
 */
double DVector::dot(const DVector& vec)
{
  int i, N;
  double dot;

  N = vec.size();
  dot = 0.0;

  for (i = 0 ; i < N ; ++i) {
    dot += this->operator[](i) * vec[i];
  }

  return dot;
}


void DVector::print() const
{
  int i;

  for (i = 0 ; i < this->size() ; ++i) {
    std::cout << this->operator[](i) << " ";
  }

  std::cout << std::endl << std::endl;
}


// Operators
/*!
 * @param [in] os the output stream in which to write the DVector.
 * @param [in] v the Dvector that we want to write.
 *
 * @return a reference to os (so that we can chain the << operators).
 *
 * @deprecated This operator was only used for debugging purposes. It is not used
 * in the final version of the code.
 */
std::ostream& operator<< (std::ostream &os, const DVector& v)
{
  int i, N;

  N = v.size();

  for (i = 0 ; i < N - 1 ; ++i) {
    os << v[i] << " ";
  }
  os << v[N-1] << std::endl;

  return os;
}


/*!
 * @param [in] u The first DVector operand.
 * @param [in] v The second DVector operand.
 *
 * @return A DVector that is the sum of u and v.
 */
DVector operator+ (const DVector& u, const DVector& v)
{
  DVector w;
  int i;

  if (u.size() != v.size()) {
    std::cout << "ERROR : DVector sizes do not match (" << u.size() << " and " << v.size() << ")" << std::endl;
    exit(EXIT_FAILURE);
  }

  w.resize(u.size());
  for (i = 0 ; i < u.size() ; ++i) {
    w[i] = u[i] + v[i];
  }

    return w;
}


/*!
 * @param [in] u The first DVector operand.
 * @param [in] v The second DVector operand.
 *
 * @return A DVector that is the difference of u and v.
 */
DVector operator- (const DVector& u, const DVector& v)
{
  DVector w;
  int i;

  if (u.size() != v.size()) {
    std::cout << "ERROR : DVector sizes do not match (" << u.size() << " and " << v.size() << ")" << std::endl;
    exit(EXIT_FAILURE);
  }

  w.resize(u.size());
  for (i = 0 ; i < u.size() ; ++i) {
    w[i] = u[i] - v[i];
  }

  return w;

}


/*!
 * @param [in] alpha The scalar to multiply the DVector with
 * @param [in] v The DVector that is multiplied.
 *
 * @return A DVector equal to \f$ \alpha v \f$.
 */
DVector operator* (double alpha, const DVector& u)
{
  DVector w;
  int i;

  w.resize(u.size());

  for (i = 0 ; i < u.size() ; ++i) {
    w[i] = alpha * u[i];
  }

  return w;
}


/*!
 * @param [in] alpha The scalar to multiply the DVector with
 * @param [in] v The DVector that is multiplied.
 *
 * @return A DVector equal to \f$ \alpha v \f$.
 */
DVector operator* (const DVector& u, double alpha)
{
  DVector w;
  int i;

  w.resize(u.size());

  for (i = 0 ; i < u.size() ; ++i) {
    w[i] = alpha * u[i];
  }

  return w;
}
