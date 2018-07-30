//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_hardening_laws/hardening_law.hpp"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(SOLID_MECHANICS_APPLICATION) JohnsonCookThermalHardeningLaw
	: public HardeningLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of JohnsonCookThermalHardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookThermalHardeningLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    JohnsonCookThermalHardeningLaw();


    /// Copy constructor.
    JohnsonCookThermalHardeningLaw(JohnsonCookThermalHardeningLaw const& rOther);

    /// Assignment operator.
    JohnsonCookThermalHardeningLaw& operator=(JohnsonCookThermalHardeningLaw const& rOther);

    /// Destructor.
    ~JohnsonCookThermalHardeningLaw();

    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this hardening law
     */
    virtual HardeningLaw::Pointer Clone() const;


    ///@}
    ///@name Operations
    ///@{

    double& CalculateHardening(double &rHardening, const Parameters& rValues);

    double& CalculateIsotropicHardening(double &rIsotropicHardening, const Parameters& rValues);

    double& CalculateKinematicHardening(double &rKinematicHardening, const Parameters& rValues);


    double& CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues);


    double& CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const Parameters& rValues);

    double& CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening, const Parameters& rValues);


    double& CalculateDeltaThermalHardening(double &rDeltaThermalHardening, const Parameters& rValues);


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    // /// Turn back information as a string.
    // std::string Info() const;

    // /// Print information about this object.
    // void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // void PrintData(std::ostream& rOStream) const;


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

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class JohnsonCookThermalHardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// // input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   JohnsonCookThermalHardeningLaw& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const JohnsonCookThermalHardeningLaw& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED  defined