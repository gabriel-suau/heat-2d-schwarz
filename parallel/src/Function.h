/*!
 * @file Function.h
 *
 * @brief Defines classes for the source terms and the boundary conditions.
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


#ifndef FUNCTION_H
#define FUNCTION_H

#include "MPIUtils.h"
#include "DataFile.h"
#include "termcolor.h"
#include "Vector.h"

/*!
 * @class Function
 *
 * @brief Handles several important functions fo the problem.
 *
 * @details This class is used to compute the initial condition, the source term, the
 * exact solution of the problem if it exists, and to handle the boundary conditions.
 */
class Function
{
protected:
  // Pointeur vers le fichier de paramètres
  DataFile* _DF; ///< Pointer to a DataFile object.

  // Variables utiles
  double _xmin; ///< xmin
  double _ymin; ///< ymin
  double _xmax; ///< xmax
  double _ymax; ///< ymax
  double _Lx; ///< Length of the domain in the x direction
  double _Ly; ///< Length of the domain in the y direction
  double _dx; ///< Space step in the x direction
  double _dy; ///< Space stepp in the y direction
  int _Nx; ///< Number of nodes in the x direction
  int _Ny; ///< Number of unknown in the y direction
  int _N;

  // Vecteur solution initiale
  DVector _Sol0; ///< Initial condition.
  // Terme source
  DVector _sourceTerm; ///< Source term.
  // Solution exacte
  DVector _exactSol; ///< Exact solution.

public:
  /*! @brief Constructs a Function object using a DataFile object. */
  Function(DataFile* DF);

  /*! @brief Default destructor. */
  virtual ~Function() = default;

  /*! @brief Initializes an already constructed Function object. */
  void Initialize();

  /*! @brief Builds the source term. */
  void buildSourceTerm(double t);
  /*! @brief Update the source term with the overlaping terms */
  void updateSourceTerm(const DVector& sol);

  /*! @brief Builds the exact solution of the selected scenario if it exists (pure virtual). */
  virtual void buildExactSolution(double t) = 0;

  /*! @brief Saves the current exact solution in a file. */
  void saveCurrentExactSolution(std::string& fileName) const;

  /*! @brief Saves the current source term in a file. */
  void saveSourceTerm(std::string& fileName) const;

  const DVector& getInitialCondition() const {return _Sol0;};
  const DVector& getSourceTerm() const {return _sourceTerm;};
  const DVector& getExactSolution() const {return _exactSol;};

  /*! @brief Evaluates the source term at point (x, y) at time t depending on the selected scenario (pure virtual). */
  virtual double f(const double x, const double y, const double t) = 0;

  /*! @brief Evaluates the Dirichlet left/right boundary condition at point (x, y) at time t depending on the selected scenario (pure virtual). */
  virtual double g(const double x, const double y, const double t) = 0;

  /*! @brief Evaluates the Dirichlet top/bottom boundary condition at point (x, y) at time t depending on the selected scenario (pure virtual). */
  virtual double h(const double x, const double y, const double t) = 0;
};

/*!
 * @class Function1
 *
 * @brief Handles the source term, the boundary conditions and the exact solution for scenario 1.
 */
class Function1 : public Function
{
public:
  // Constructeur
  Function1(DataFile* DF);

  // Fonctions
  double f(const double x, const double y, const double t);
  double g(const double x, const double y, const double t);
  double h(const double x, const double y, const double t);

  // Exact solution
  void buildExactSolution(double t);
};

/*!
 * @class Function2
 *
 * @brief Handles the source term, the boundary conditions and the exact solution for scenario 2.
 */
class Function2 : public Function
{
public:
  // Constructeur
  Function2(DataFile* DF);

  // Fonctions
  double f(const double x, const double y, const double t);
  double g(const double x, const double y, const double t);
  double h(const double x, const double y, const double t);

  // Exact solution
  void buildExactSolution(double t);
};

/*!
 * @class Function3
 *
 * @brief Handles the source term and the boundary conditions for scenario 3.
 */
class Function3 : public Function
{
public:
  // Constructeur
  Function3(DataFile* DF);

  // Fonctions
  double f(const double x, const double y, const double t);
  double g(const double x, const double y, const double t);
  double h(const double x, const double y, const double t);

  // Exact solution
  void buildExactSolution(double t);
};

#endif // FUNCTION_H
