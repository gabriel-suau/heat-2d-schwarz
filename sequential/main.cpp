#include "termcolor.h"
#include "Vector.h"
#include "DataFile.h"
#include "Function.h"
#include "Laplacian.h"
#include "TimeScheme.h"

#include <iostream>
#include <chrono>

int main(int argc, char** argv)
{
  //------------------------------------------------//
  //---------------------Checks---------------------//
  //------------------------------------------------//
  if (argc < 2)
    {
      std::cout << termcolor::red << "ERROR::MAIN : Please, enter the name of your data file." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
  else if (argc > 2)
    {
      std::cout << termcolor::yellow << "WARNING::MAIN : Too many arguments : only 1 was expected but " << argc - 1 << " were given..."<< std::endl;
      std::cout << termcolor::reset;
    }
  
  //--------------------------------------------------------//
  //---------------------Beginning logs---------------------//
  //--------------------------------------------------------//
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Solving 2D heat equation for you !" << std::endl;
  std::cout << "====================================================================================================" << std::endl << std::endl;

  
  //--------------------------------------------------------//
  //---------------------Parameters (Data) File-------------//
  //--------------------------------------------------------//
  DataFile* DF = new DataFile(argv[1]);
  DF->readDataFile();
  DF->printData();

  
  //-------------------------------------------------------------//
  //---------------------IC, BC, Source term---------------------//
  //-------------------------------------------------------------//
  Function* function;
  if (DF->getScenario() == 1)
    {
      function = new Function1(DF);
    }
  else if (DF->getScenario() == 2)
    {
      function = new Function2(DF);
    }
  else if (DF->getScenario() == 3)
    {
      function = new Function3(DF);
    }
  else
    {
      std::cout << termcolor::red << "ERROR::FUNCTION : Scenario not implemented !" << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
  function->Initialize();
  

  //------------------------------------------------------------//
  //---------------------Discrete laplacian---------------------//
  //------------------------------------------------------------//
  Laplacian* laplacian = new Laplacian(DF, function);
  laplacian->Initialize();
  
  
  //-----------------------------------------------------//
  //---------------------Time Scheme---------------------//
  //-----------------------------------------------------//
  TimeScheme* TS;
  if (DF->getTimeScheme() == "ExplicitEuler")
    {
      TS = new ExplicitEuler(DF, function, laplacian);
    }
  else if (DF->getTimeScheme() == "ImplicitEuler")
    {
      TS = new ImplicitEuler(DF, function, laplacian);
    }
  else
    {
      std::cout << termcolor::red << "ERROR::TIMESCHEME : Case not implemented." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
  

  //-------------------------------------------------//
  //---------------------Solving---------------------//
  //-------------------------------------------------//
  TS->solve();
  
  
  //---------------------------------------------------------//
  //---------------------Free the memory---------------------//
  //---------------------------------------------------------//
  delete DF;
  delete function;
  delete laplacian;
  delete TS;

  
  //-----------------------------------------------------//
  //---------------------Ending logs---------------------//
  //-----------------------------------------------------//
  std::cout << "====================================================================================================" << std::endl;
  std::cout << termcolor::green << "SUCCESS : Successfully solved the 2D heat equation for you !" << std::endl;
  std::cout << termcolor::reset << "Let me terminate myself now..." << std::endl;
  std::cout << "====================================================================================================" << std::endl << std::endl;

    
  // End
  return EXIT_SUCCESS;
}
