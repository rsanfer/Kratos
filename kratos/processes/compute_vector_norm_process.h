//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined(KRATOS_COMPUTE_VECTOR_NORM_PROCESS_H_INCLUDED )
#define  KRATOS_COMPUTE_VECTOR_NORM_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/enums.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Compute the norm of a vector
/** This process computes the norm of a cartain vector variable stored in the nodes
*/
template< HistoricalValues THist>
class KRATOS_API(KRATOS_CORE) ComputeVectorNormProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeVectorNormProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeVectorNormProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeVectorNormProcess(
        ModelPart& rModelPart,
        Variable<array_1d<double,3> >& rOriginVariable,
        Variable<double>& rNormVariable);

    /// Destructor.
    virtual ~ComputeVectorNormProcess() override
    {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    /**
     * Execute method is used to execute the Process algorithms.
     * In this process the gradient of a scalar variable will be computed
     */
    void Execute() override;


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ComputeVectorNormProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeVectorNormProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Variable<array_1d<double, 3> >& mrOriginVariable;
    Variable<double>& mrNormVariable;

    ///@}
    ///@name Private Operators
    ///@{

    void ComputeNorm();

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeVectorNormProcess& operator=(ComputeVectorNormProcess const& rOther);

    /// Copy constructor.
    ComputeVectorNormProcess(ComputeVectorNormProcess const& rOther);


  ///@}

}; // Class ComputeVectorNormProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
// 			    ComputeVectorNormProcess& rThis);

/// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
// 			    const ComputeVectorNormProcess& rThis)
// {
//   rThis.PrintInfo(rOStream);
//   rOStream << std::endl;
//   rThis.PrintData(rOStream);
//
//   return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COMPUTE_VECTOR_NORM_PROCESS_H_INCLUDED  defined
