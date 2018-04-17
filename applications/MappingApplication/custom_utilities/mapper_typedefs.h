//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes


// Project includes
#include "spaces/ublas_space.h"
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "trilinos_space.h"
#include "Epetra_FEVector.h"
// #include "Epetra_FECrsMatrix.h" // included in "trilinos_space"
#endif

namespace Kratos
{
    namespace MapperDefinitions
    {
        typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

        typedef UblasSpace<double, CompressedMatrix, Vector> UblasSparseSpaceType;

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
#endif
    }

}  // namespace Kratos.