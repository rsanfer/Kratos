//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotoamyor
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/variable_utils.h"
#include "processes/compute_vector_norm_process.h"

namespace Kratos
{

template< HistoricalValues THist >
ComputeVectorNormProcess<THist>::ComputeVectorNormProcess(
    ModelPart& rModelPart,
    Variable<array_1d<double,3> >& rOriginVariable,
    Variable<double>& rNormVariable)
    : mrModelPart(rModelPart),
      mrOriginVariable(rOriginVariable),
      mrNormVariable(rNormVariable)
{
    KRATOS_TRY

    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rNormVariable, mrModelPart.Nodes());

    KRATOS_CATCH("")
}



}
