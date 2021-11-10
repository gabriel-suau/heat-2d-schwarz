
/*!
 * @file DataFile.cpp
 *
 * @brief Define a class representing a Data file (contains all the parameters of the simulation).
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


#include "DataFile.h"
#include "termcolor.h"
#include "MPIUtils.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <regex>
#include <mpi.h>


DataFile::DataFile(const std::string& fileName):
  _fileName(fileName)
{
}


/*!
 * @details Removes everything after a possible #, replaces tabulation by spaces,
 * replaces multiple spaces by one space, removes any leading space.
 *
 * @param [in] line The line to clean.
 *
 * @return A std::string corresponding to the cleaned line.
 */
std::string DataFile::cleanLine(std::string &line)
{
  std::string res = line;

  // Remove everything after a possible #
  res = regex_replace(res, std::regex("#.*$"), std::string(""));
  // Replace tabulation(s) by space(s)
  res = regex_replace(res, std::regex("\t"), std::string(" "), std::regex_constants::match_any);
  // Replace multiple spaces by 1 space
  res = regex_replace(res, std::regex("\\s+"), std::string(" "), std::regex_constants::match_any);
  // Remove any leading spaces
  res = regex_replace(res, std::regex("^ *"), std::string(""));

  return res;
}


/*!
 * @details Main method of this class. It reads the data file fileName, and assign the values to the
 * parameters. It also performs checks to ensure the values entered by the user in the data file
 * are correct.
 */
void DataFile::readDataFile()
{
  // Open the data file
  std::ifstream dataFile(_fileName.data());
  if (!dataFile.is_open())
    {
      if (MPI_Rank == 0)
        {
          std::cout << termcolor::red << "ERROR::DATAFILE : Unable to open file " << _fileName << std::endl;
          std::cout << termcolor::reset;
        }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
#if VERBOSITY>0
  else
    {
      if (MPI_Rank == 0)
        {
          std::cout << "====================================================================================================" << std::endl;
          std::cout << "Reading data file " << _fileName << std::endl;
        }
    }
#endif

  // Pour stocker chaque ligne
  std::string line;
  // Run through the dataFile to find the parameters
  while (getline(dataFile, line))
    {
      // Clean line
      std::string proper_line(cleanLine(line));
      if (proper_line.find("Scenario") != std::string::npos)
        {
          dataFile >> _scenario;
        }
      if (proper_line.find("ResultsDir") != std::string::npos)
        {
          dataFile >> _resultsDir;
        }
      if (proper_line.find("SaveFinalResultOnly") != std::string::npos)
        {
          dataFile >> _isSaveFinalResultOnly;
        }
      if (proper_line.find("SaveFrequency") != std::string::npos)
        {
          dataFile >> _saveFrequency;
        }
      if (proper_line.find("ErrorAndCPUTimeDir") != std::string::npos)
        {
          dataFile >> _errorAndCPUTimeDir;
        }
      if (proper_line.find("TimeScheme") != std::string::npos)
        {
          dataFile >> _timeScheme;
        }
      if (proper_line.find("InitialTime") != std::string::npos)
        {
          dataFile >> _initialTime;
        }
      if (proper_line.find("FinalTime") != std::string::npos)
        {
          dataFile >> _finalTime;
        }
      if (proper_line.find("TimeStep") != std::string::npos)
        {
          dataFile >> _timeStep;
        }
      if (proper_line.find("CFL") != std::string::npos)
        {
          dataFile >> _CFL;
        }
      if (proper_line.find("nSubdomains") != std::string::npos)
        {
          dataFile >> _nSubdomains;
        }
      if (proper_line.find("nOverlap") != std::string::npos)
        {
          dataFile >> _nOverlap;
        }
      if (proper_line.find("xmin") != std::string::npos)
        {
          dataFile >> _xmin;
        }
      if (proper_line.find("xmax") != std::string::npos)
        {
          dataFile >> _xmax;
        }
      if (proper_line.find("ymin") != std::string::npos)
        {
          dataFile >> _ymin;
        }
      if (proper_line.find("ymax") != std::string::npos)
        {
          dataFile >> _ymax;
        }
      if (proper_line.find("Nx") != std::string::npos)
        {
          dataFile >> _Nx;
        }
      if (proper_line.find("Ny") != std::string::npos)
        {
          dataFile >> _Ny;
        }
      if (proper_line.find("MaxIterations") != std::string::npos)
        {
          dataFile >> _maxIterations;
        }
      if (proper_line.find("Tolerance") != std::string::npos)
        {
          dataFile >> _tolerance;
        }
      if (proper_line.find("IsSaveResidual") != std::string::npos)
        {
          dataFile >> _isSaveResidual;
        }
      if (proper_line.find("ResidualFile") != std::string::npos)
        {
          dataFile >> _resFile;
        }
      if (proper_line.find("DiffusionCoefficient") != std::string::npos)
        {
          dataFile >> _diffCoeff;
        }
    }

  // Calcul des pas d'espace
  _Lx = _xmax - _xmin;
  _Ly = _ymax - _ymin;
  _dx = _Lx / (_Nx + 1);
  _dy = _Ly / (_Ny + 1);

  // Calcul du pas de temps pour Euler Explicite
  if (_timeScheme == "ExplicitEuler")
    {
#if VERBOSITY>0
      if (MPI_Rank == 0)
        {
          std::cout << termcolor::yellow << "Adjusting the time step to fit the CFL condition" << std::endl;
          std::cout << termcolor::reset;
        }
#endif
      _timeStep = _CFL * (pow(_dx,2) * pow(_dy,2)) / (2. * _diffCoeff * (pow(_dx,2) + pow(_dy,2)));
    }

  // Calcul du nombre d'iterations en temps et ajustement du pas de temps
#if VERBOSITY>0
  if (MPI_Rank == 0)
    {
      std::cout << termcolor::yellow << "Adjusting the time step to land exactly on the final time" << std::endl;
      std::cout << termcolor::reset;
    }
#endif

  int nbIterations(int(ceil((_finalTime - _initialTime)/_timeStep)));
  _timeStep = (_finalTime - _initialTime)/nbIterations;

#if VERBOSITY>0
  if (MPI_Rank == 0)
    std::cout << "The new time step is dt = " << _timeStep << std::endl;
#endif

  // Création et nettoyage du dossier de résultats
  if (MPI_Rank == 0)
    {

#if VERBOSITY>0
      std::cout << "Creating the results directory..." << std::endl;
      std::cout << "Creating the error and cputime directory..." << std::endl;
#endif

      // Results directory
      system(("mkdir -p ./" + _resultsDir).c_str());
      system(("rm -f ./" + _resultsDir + "/solution*").c_str());
      system(("rm -f ./" + _resultsDir + "/" + _resFile).c_str());
      system(("cp -r ./" + _fileName + " ./" + _resultsDir + "/params.txt").c_str());
      // Error and cputime directory
      system(("mkdir -p ./" + _errorAndCPUTimeDir).c_str());

#if VERBOSITY>0
      // Logs
      std::cout << termcolor::green << "SUCCESS::DATAFILE : Results directory created successfully !" << std::endl;
      std::cout << termcolor::reset;

      // Logs de succès
      std::cout << termcolor::green << "SUCCESS::DATAFILE : File read successfully" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
#endif
    }
}


// Affiche les paramètres sur le terminal
void DataFile::printData() const
{
#if VERBOSITY>0
  if (MPI_Rank == 0)
    {
      std::cout << "====================================================================================================" << std::endl;
      std::cout << "Printing parameters of " << _fileName << std::endl;
      std::cout << "Scenario              = " << _scenario << std::endl;
      std::cout << "Diffusion Coefficient = " << _diffCoeff << std::endl;
      std::cout << "Spatial parameters : " << std::endl;
      std::cout << "    |xmin             = " << _xmin << std::endl;
      std::cout << "    |xmax             = " << _xmax << std::endl;
      std::cout << "    |ymin             = " << _ymin << std::endl;
      std::cout << "    |ymax             = " << _ymax << std::endl;
      std::cout << "    |Lx               = " << _Lx << std::endl;
      std::cout << "    |Ly               = " << _Ly << std::endl;
      std::cout << "    |Nx               = " << _Nx << std::endl;
      std::cout << "    |Ny               = " << _Ny << std::endl;
      std::cout << "    |dx               = " << _dx << std::endl;
      std::cout << "    |dy               = " << _dy << std::endl;
      std::cout << "    |Nb domains       = " << _nSubdomains << std::endl;
      std::cout << "    |Overlap size     = " << _nOverlap << std::endl;
      std::cout << "Time Scheme           = " << _timeScheme << std::endl;
      std::cout << "Initial time          = " << _initialTime << std::endl;
      std::cout << "Final time            = " << _finalTime << std::endl;
      if (_timeScheme == "ExplicitEuler")
        {
          std::cout << "CFL number            = " << _CFL << std::endl;
        }
      std::cout << "Time step             = " << _timeStep << std::endl;
      if (_timeScheme == "ImplicitEuler")
        {
          std::cout << "Linear solver         : Conjugate Gradient" << std::endl;
          std::cout << "    |Max Iterations   = " << _maxIterations << std::endl;
          std::cout << "    |Tolerance        = " << _tolerance << std::endl;
          std::cout << "    |Save residual ?  = " << _isSaveResidual << std::endl;
          if (_isSaveResidual)
            std::cout << "    |Residual File    = " << _resFile << std::endl;
        }
      std::cout << "Results directory     = " << _resultsDir << std::endl;
      std::cout << "Save final time only? = " << _isSaveFinalResultOnly << std::endl;
      if (!_isSaveFinalResultOnly)
        std::cout << "Save Frequency        = " << _saveFrequency << std::endl;
      std::cout << "====================================================================================================" << std::endl << std::endl;
    }
#endif
}
