#ifndef TIME_SCHEME_H
#define TIME_SCHEME_H

#include "DataFile.h"
#include "Function.h"
#include "Laplacian.h"
#include "Vector.h"
#include <string>

class TimeScheme
{
protected:
  // Pointeur vers les trucs importants
  DataFile* _DF;
  Function* _function;
  Laplacian* _laplacian;

  // Solution
  DVector _Sol;
  
  // Paramètres de temps
  double _timeStep;
  double _initialTime;
  double _finalTime;
  double _currentTime;

  // Sauvegarde des résultats
  std::string _resultsDir;
  std::string _resFileName;
  
public:
  // Constructeurs
  TimeScheme();
  TimeScheme(DataFile* DF, Function* function, Laplacian* laplacian);

  // Initialiseur
  void Initialize(DataFile* DF, Function* function, Laplacian* laplacian);
  // Destructeur
  virtual ~TimeScheme() = default;

  // Getters
  const DVector& getSolution() const {return _Sol;};
  double getTimeStep() const {return _timeStep;};
  double getInitialTime() const {return _initialTime;};
  double getFinalTime() const {return _finalTime;};
  double getCurrentTime() const {return _currentTime;};
  
  // Solve and save solution
  virtual void oneStep() = 0;
  void saveCurrentSolution(std::string& fileName) const;
  void solve();

  // Compute the L2 error norm
  double computeCurrentError();
};

class ExplicitEuler: public TimeScheme
{
public:
  // Constructeurs
  ExplicitEuler();
  ExplicitEuler(DataFile* DF, Function* function, Laplacian* laplacian);

  // Initialiseur
  void Initialize(DataFile* DF, Function* function, Laplacian* laplacian);

  // One time step
  void oneStep();
};


class ImplicitEuler: public TimeScheme
{
public:
  // Constructeurs
  ImplicitEuler();
  ImplicitEuler(DataFile* DF, Function* function, Laplacian* laplacian);

  // Initialiseur
  void Initialize(DataFile* DF, Function* function, Laplacian* laplacian);

  // One time step
  void oneStep();
};

#endif // TIME_SCHEME_H
