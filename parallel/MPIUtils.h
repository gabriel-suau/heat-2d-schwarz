/*!
 * @file MPIUtils.h
 *
 * @brief Defines the global MPI variables and a function to compute the load allocated to each MPI proc.
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

#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <mpi.h>

// Global MPI variables
extern int MPI_Size; ///< Total number of of MPI procs.
extern int MPI_Rank; ///< Rank of this MPI proc.
extern int kBegin; ///< Index of the first unknown allocated to this MPI proc.
extern int kEnd; ///< Index of the last unknown allocated to this MPI proc.
extern int localSize; ///< Number of unknown allocated to this MPI proc.
extern int rowBegin; ///< Index of the first row of the domain allocated to this MPI proc.
extern int rowEnd; ///< Index of the last row of the domain allocated to this MPI proc.
extern int nbDomainRows; ///< Number of rows of the domain allocated to this MPI proc.
extern MPI_Status status; ///< Status of the communicator.


/*!
 * @brief Computes the load allocated to an MPI proc.
 *
 * @details Compute the load allocated to an MPI proc knowing the total load, the number of MPI procs,
 * and the rank of this MPI proc, with a maximum load imbalance of 1.
 *
 * @param[in] N Total load to allocate.
 * @param[in] Np Number of MPI procs.
 * @param[in] me Rank of this MPI proc.
 * @param[out] iBegin Index of the first element allocated to this MPI proc.
 * @param[out] iEnd Index of the last element allocated to this MPI proc.
 */
inline void charge(int N, int Np, int me, int nOverlap, int* iBegin, int* iEnd)
{
  // Division entière
  int chargeMinParProc(N/Np);
  // Reste à répartir
  int reste(N%Np);
  // Gestion du recouvrement
  // nOverlap pair -->
  // nOverlap impair -->
  int add_up, add_down;
  add_up = nOverlap / 2;
  add_down = add_up + nOverlap % 2;

  if (me < reste)
    {
      *iBegin = me * (chargeMinParProc + 1);
      *iEnd = *iBegin + chargeMinParProc;
    }
  else
    {
      *iBegin = reste + me * chargeMinParProc;
      *iEnd = *iBegin + chargeMinParProc - 1;
    }

  if (me + 1 < MPI_Size) {
    *iEnd += add_up;
  }
  if (me > 0) {
    *iBegin -= add_down;
  }

}

#endif // MPI_UTILS_H
