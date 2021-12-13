/*!
 * @file DataFile.h
 *
 * @brief DataFile class to read the simulation parameters.
 *
 * This file contains a DataFile class that reads the parameters file
 * and contains all the parameters of the simulation.
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


#ifndef DATA_FILE_H
#define DATA_FILE_H

#include "MPIUtils.h"
#include <string>

/*!
 * @class DataFile
 *
 * @brief Represents a data file.
 *
 * @details A DataFile object contains all the parameters of the simulation. The parameters are read
 * in a formatted data file.
*/
class DataFile
{
private:
  std::string _fileName; ///< Name of the data file.

  // Solution saving
  std::string _resultsDir; ///< Directory in which the results are written.
  bool _isSaveFinalResultOnly; ///< Boolean to check wether only the final result must be saved.
  int _saveFrequency; ///< Number of iterations between each saving of the solution (ignored if _isSaveFinalResultOnly is set to true).
  std::string _errorAndCPUTimeDir; ///< Directory in which the CPU time and the error are written.

  // Scenario
  int _scenario; ///< Scenario that is simulated (1, 2 or 3).

  // Time parameters
  std::string _timeScheme; ///< Time integrator (ExplicitEuler or ImplicitEuler).
  double _initialTime; ///< Initial time.
  double _finalTime; ///< Final time.
  double _timeStep; ///< Time step.
  double _CFL; ///< CFL number. It is used with the Explicit Euler scheme to ensure its stability.

  // Spatial parameters
  double _xmin, _xmax, _ymin, _ymax; ///< Boundaries of the domain.
  double _Lx, _Ly; ///< Lengths of the domain in the x and y directions.
  int _Nx, _Ny; ///< Number of nodes (or unknowns) in the x and y directions.
  double _dx, _dy; ///< Space steps in the x and y directions.

  // Schwarz parameters
  int _nOverlap; ///< Number of subdomains, number of over lap lines
  int _schwarzMaxIterations;
  double _schwarzTolerance;

  // Linear solver parameters
  std::string _linearSolver;
  int _maxIterations; ///< Maximum number of iterations for the Conjugate Gradient.
  double _tolerance; ///< Tolerance for the Conjugate Gradient.

  // Diffusion coefficient
  double _diffCoeff; ///< Diffusion coefficient

public:
  /*! @brief Construct a DataFile object using the data file name. */
  DataFile(const std::string& fileName);

  /*! @brief Default destructor. */
  ~DataFile() = default;

  /*! @brief Reads the data file and sets the values of the parameters. */
  void readDataFile();

  // Getters

  // DataFile name
  const std::string& getFileName() const {return _fileName;}

  // Solution saving
  const std::string& getResultsDirectory() const {return _resultsDir;}
  bool isSaveFinalResultOnly() const {return _isSaveFinalResultOnly;}
  int getSaveFrequency() const {return _saveFrequency;}
  const std::string& getErrorAndCPUTimeDir() const {return _errorAndCPUTimeDir;}

  // Scenario
  int getScenario() const {return _scenario;}

  // Time parameters
  const std::string& getTimeScheme() const {return _timeScheme;}
  double getInitialTime() const {return _initialTime;}
  double getFinalTime() const {return _finalTime;}
  double getTimeStep() const {return _timeStep;}
  double getCFL() const {return _CFL;}

  // Mesh parameters
  double getxMin() const {return _xmin;}
  double getxMax() const {return _xmax;}
  double getyMin() const {return _ymin;}
  double getyMax() const {return _ymax;}
  double getLx() const {return _Lx;}
  double getLy() const {return _Ly;}
  int getNx() const {return _Nx;}
  int getNy() const {return _Ny;}
  double getDx() const {return _dx;}
  double getDy() const {return _dy;}

  // Schwarz parameters
  int getnOverlap() const {return _nOverlap;}
  int getSchwarzMaxIterations() const {return _schwarzMaxIterations;}
  double getSchwarzTolerance() const {return _schwarzTolerance;}

  // Linear solver parameters
  std::string getLinearSolver() const {return _linearSolver;}
  int getMaxIterations() const {return _maxIterations;}
  double getTolerance() const {return _tolerance;}

  // Diffusion coefficient
  double getDiffCoeff() const {return _diffCoeff;}

  /*! @brief Prints the values of the parameters. */
  void printData() const;

protected:
  /*! @brief Cleans a line so that it can be read. */
  std::string cleanLine(std::string &line);
};

#endif // DATA_FILE_H
