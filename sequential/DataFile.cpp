#include "DataFile.h"
#include "termcolor.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <regex>

DataFile::DataFile()
{
}

DataFile::DataFile(const std::string &fileName) : _fileName(fileName)
{
}

void DataFile::Initialize(const std::string &fileName)
{
  _fileName = fileName;
}

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

void DataFile::readDataFile()
{
  // Open the data file
  std::ifstream dataFile(_fileName.data());
  if (!dataFile.is_open())
  {
    std::cout << termcolor::red << "ERROR::DATAFILE : Unable to open file " << _fileName << std::endl;
    std::cout << termcolor::reset;
    exit(-1);
  }
  else
  {
    std::cout << "====================================================================================================" << std::endl;
    std::cout << "Reading data file " << _fileName << std::endl;
  }
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
    if (proper_line.find("DiffusionCoefficient") != std::string::npos)
    {
      dataFile >> _diffCoeff;
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
    if (proper_line.find("ResultsDir") != std::string::npos)
    {
      dataFile >> _resultsDir;
    }
    if (proper_line.find("SaveFrequency") != std::string::npos)
    {
      dataFile >> _saveFrequency;
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
    std::cout << termcolor::yellow << "Adjusting the time step to fit the CFL condition" << std::endl;
    std::cout << termcolor::reset;
    _timeStep = _CFL * (pow(_dx, 2) * pow(_dy, 2)) / (2. * _diffCoeff * (pow(_dx, 2) + pow(_dy, 2)));
  }

  // Calcul du nombre d'iterations en temps et ajustement du pas de temps
  std::cout << termcolor::yellow << "Adjusting the time step to land exactly on the final time" << std::endl;
  std::cout << termcolor::reset;
  int nbIterations(int(ceil((_finalTime - _initialTime) / _timeStep)));
  _timeStep = (_finalTime - _initialTime) / nbIterations;
  std::cout << "The new time step is dt = " << _timeStep << std::endl;

  // Création et nettoyage du dossier de résultats
  std::cout << "Creating the results directory..." << std::endl;
  system(("mkdir -p ./" + _resultsDir).c_str());
  system(("rm -f ./" + _resultsDir + "/solution*").c_str());
  system(("rm -f ./" + _resultsDir + "/" + _resFile).c_str());
  system(("cp -r ./" + _fileName + " ./" + _resultsDir + "/params.txt").c_str());

  // Logs
  std::cout << termcolor::green << "SUCCESS::DATAFILE : Results directory created successfully !" << std::endl;
  std::cout << termcolor::reset;

  // Logs de succès
  std::cout << termcolor::green << "SUCCESS::DATAFILE : File read successfully" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl
            << std::endl;
}

// Affiche les paramètres sur le terminal
void DataFile::printData() const
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
    {
      std::cout << "    |Residual File    = " << _resFile << std::endl;
    }
  }
  std::cout << "Results directory     = " << _resultsDir << std::endl;
  std::cout << "Save Frequency        = " << _saveFrequency << std::endl;
  std::cout << "====================================================================================================" << std::endl
            << std::endl;
}
