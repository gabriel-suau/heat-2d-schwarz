#ifndef DATA_FILE_H
#define DATA_FILE_H

#include <string>

class DataFile
{
private:
  // Name of the data file
  std::string _fileName;

  // Directory in which the results are saved
  std::string _resultsDir;

  // Initial condition
  int _scenario;

  // Time parameters
  std::string _timeScheme;
  double _initialTime;
  double _finalTime;
  double _timeStep;
  double _CFL;

  // Spatial parameters
  double _xmin, _xmax, _ymin, _ymax;
  double _Lx, _Ly;
  int _Nx, _Ny;
  double _dx, _dy;

  int _nSubdomains;
  int _nOverlapLine;

  // Linear solver parameters
  int _maxIterations;
  double _tolerance;
  bool _isSaveResidual;
  std::string _resFile;

  // Diffusion coefficient
  double _diffCoeff;

  // Frequency at which results are saved
  int _saveFrequency;

public:
  DataFile();
  DataFile(const std::string &fileName);

  ~DataFile() = default;

  void Initialize(const std::string &fileName);

  void readDataFile();

  std::string cleanLine(std::string &line);

  // Getters
  const std::string &getFileName() const { return _fileName; };
  const std::string &getResultsDirectory() const { return _resultsDir; };
  int getScenario() const { return _scenario; };
  const std::string &getTimeScheme() const { return _timeScheme; };
  double getInitialTime() const { return _initialTime; };
  double getFinalTime() const { return _finalTime; };
  double getTimeStep() const { return _timeStep; };
  double getCFL() const { return _CFL; };
  double getxMin() const { return _xmin; };
  double getxMax() const { return _xmax; }
  double getyMin() const { return _ymin; };
  double getyMax() const { return _ymax; };
  double getLx() const { return _Lx; };
  double getLy() const { return _Ly; };
  int getNx() const { return _Nx; };
  int getNy() const { return _Ny; };
  double getDx() const { return _dx; }
  double getDy() const { return _dy; };
  int getMaxIterations() const { return _maxIterations; };
  double getTolerance() const { return _tolerance; };
  bool isSaveResidual() const { return _isSaveResidual; }
  std::string getResFile() const { return _resFile; };
  double getDiffCoeff() const { return _diffCoeff; };
  int getSaveFrequency() const { return _saveFrequency; };
  int getnSubsomains() const { return _saveFrequency; };
  int getnOverlapLines() const { return _saveFrequency; };

  // Print the parameters
  void printData() const;
};

#endif // DATA_FILE_H
