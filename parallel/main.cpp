/*!
 * @file main.cpp
 *
 * @brief Main file of the program.
 *
 * @details Initializes the MPI Environment, and calls the functions to solve the 2D heat equation.
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

#include "termcolor.h"
#include "MPIUtils.h"
#include "Vector.h"
#include "DataFile.h"
#include "Function.h"
#include "Laplacian.h"
#include "TimeScheme.h"

#include <mpi.h>
#include <iostream>

// Global MPI variables
int MPI_Size, MPI_Rank, kBegin, kEnd, rowBegin, rowEnd, localSize, nbDomainRows;
MPI_Status status;

// Main function
int main(int argc, char** argv)
{
  //------------------------------------------------------//
  //---------------------Init MPI Env---------------------//
  //------------------------------------------------------//
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_Size);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);
  

  //------------------------------------------------//
  //---------------------Checks---------------------//
  //------------------------------------------------//
  if (argc < 2)
    {
      if (MPI_Rank == 0)
        {
          std::cout << termcolor::red << "ERROR::MAIN : Please, enter the name of your data file." << std::endl;
          std::cout << termcolor::reset;
        }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  else if (argc > 2)
    {
      if (MPI_Rank == 0)
        {
          std::cout << termcolor::yellow << "WARNING::MAIN : Too many arguments : only 1 was expected but " << argc - 1 << " were given..."<< std::endl;
          std::cout << termcolor::reset; 
        }
    }


  //--------------------------------------------------------//
  //---------------------Beginning logs---------------------//
  //--------------------------------------------------------//
  if (MPI_Rank == 0)
    {
      std::cout << "====================================================================================================" << std::endl;
      std::cout << "Solving 2D heat equation for you !" << std::endl;
      std::cout << "====================================================================================================" << std::endl << std::endl; 
    }


  //--------------------------------------------------------//
  //---------------------Parameters (Data) File-------------//
  //--------------------------------------------------------//
  DataFile* DF = new DataFile(argv[1]);
  DF->readDataFile();
#if VERBOSITY>0
  DF->printData();
#endif


  //----------------------------------------------------------------//
  //---------------------Compute the load for each proc-------------//
  //----------------------------------------------------------------//
  // Compute the number of rows allocated to the proc (must be between 1 and Ny)
  charge(DF->getNx(), MPI_Size, MPI_Rank, 5, &rowBegin, &rowEnd);
  nbDomainRows = rowEnd - rowBegin + 1;
  // Local size allocated to the proc (must be between Nx and Nx*Ny)
  kBegin = DF->getNx() * rowBegin;
  kEnd = DF->getNx() * (rowEnd + 1) - 1;
  localSize = DF->getNx() * nbDomainRows;

  std::cout << MPI_Rank << " " << rowBegin << " " << rowEnd << std::endl;

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
      if (MPI_Rank == 0)
        {
          std::cout << termcolor::red << "ERROR::FUNCTION : Scenario not implemented !" << std::endl;
          std::cout << termcolor::reset; 
        }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      exit(1);
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
      if (MPI_Rank == 0)
        {
          std::cout << termcolor::red << "ERROR::TIMESCHEME : Case not implemented." << std::endl;
          std::cout << termcolor::reset; 
        }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      exit(1);
    }
  

  //-------------------------------------------------//
  //---------------------Solving---------------------//
  //-------------------------------------------------//
  // TS->solve();


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
  if (MPI_Rank == 0)
    {
      std::cout << "====================================================================================================" << std::endl;
      std::cout << termcolor::green << "SUCCESS : Successfully solved the 2D heat equation for you !" << std::endl;
      std::cout << termcolor::reset << "Let me terminate myself now..." << std::endl;
      std::cout << "====================================================================================================" << std::endl << std::endl; 
    }


  //-------------------------------------------------------//
  //---------------------Close MPI Env---------------------//
  //-------------------------------------------------------//  
  MPI_Finalize();
  
  // End
  return EXIT_SUCCESS;
}
