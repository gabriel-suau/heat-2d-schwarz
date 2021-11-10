/*!
 * @file TimeScheme.h
 *
 * @brief Defines classes for time integration.
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


#ifndef TIME_SCHEME_H
#define TIME_SCHEME_H

#include "MPIUtils.h"
#include "DataFile.h"
#include "Function.h"
#include "Laplacian.h"
#include "Vector.h"
#include <string>

/*!
 * @class TimeScheme
 *
 * @brief Represents a time integrator.
 *
 * @details A TimeScheme object contains everything that is needed to perform the time integration of the solution.
*/
class TimeScheme
{
protected:
  // Pointeur vers les trucs importants
  DataFile* _DF; ///< Pointer to a DataFile object.
  Function* _function; ///< Pointer to a Function object.
  Laplacian* _laplacian; ///< Pointer to a Laplacian object.

  // Solution
  DVector _Sol; ///< Solution vector. (updated at each time step)

  // Paramètres de temps
  double _timeStep; ///< Time step of the simulation.
  double _initialTime; ///< Initial time.
  double _finalTime; ///< Final time.
  double _currentTime; ///< Current time.

  // Sauvegarde des résultats
  std::string _resultsDir; ///< Directory in which the results are saved.

public:
  /*! @brief Constructs a TimeScheme object using a DataFile, a Function and a Laplacian objects. */
  TimeScheme(DataFile* DF, Function* function, Laplacian* laplacian);

  /*! @brief Default destructor. */
  virtual ~TimeScheme() = default;

  // Getters
  const DVector& getSolution() const {return _Sol;};
  double getTimeStep() const {return _timeStep;};
  double getInitialTime() const {return _initialTime;};
  double getFinalTime() const {return _finalTime;};
  double getCurrentTime() const {return _currentTime;};

  /*! @brief Performs one step of the selected time integration method and updates the solution (pure virtual). */
  virtual void oneStep() = 0;

  /*! @brief Saves the current solution in a file. */
  void saveCurrentSolution(std::string& fileName) const;

  /*! @brief Solves the problem. */
  void solve();

  /*! @brief Computes the current L2 error using the current solution and the current exact solution (if it exists). */
  double computeCurrentL2Error();
  /*! @brief Computes the current L1 error using the current solution and the current exact solution (if it exists). */
  double computeCurrentL1Error();
};

/*!
 * @class ExplicitEuler
 *
 * @brief Represents the Explicit Euler time scheme.
 */
class ExplicitEuler: public TimeScheme
{
public:
  ExplicitEuler(DataFile* DF, Function* function, Laplacian* laplacian);
  void oneStep();
};

/*!
 * @class ImplicitEuler
 *
 * @brief Represents the Implicit Euler time scheme.
 */
class ImplicitEuler: public TimeScheme
{
public:
  ImplicitEuler(DataFile* DF, Function* function, Laplacian* laplacian);
  void oneStep();
};

#endif // TIME_SCHEME_H
