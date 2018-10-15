// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//    Co-authors:    
//

#if !defined(KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CONSTITUTIVE_LAWS_H_INCLUDED)
#define  KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CONSTITUTIVE_LAWS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

// Advanced Constitutive laws
#include "custom_constitutive/small_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/finite_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/generic_finite_strain_isotropic_plasticity.h"
#include "custom_constitutive/generic_small_strain_isotropic_damage.h"
#include "custom_constitutive/small_strain_isotropic_damage_factory.h"
#include "custom_constitutive/viscous_generalized_kelvin.h"
#include "custom_constitutive/viscous_generalized_maxwell.h"
#include "custom_constitutive/generic_small_strain_viscoplasticity_3d.h"
#include "custom_constitutive/generic_small_strain_d_plus_d_minus_damage.h"

// HyperElastic Laws
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"

// Integrators
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_finite_strain_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/d+d-constitutive_law_integrators/generic_compression_constitutive_law_integrator.h"
#include "custom_constitutive/constitutive_laws_integrators/d+d-constitutive_law_integrators/generic_tension_constitutive_law_integrator.h"

/* Small strain */ // TODO: Move to independent folder and rename
// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

/* Finite strain */
// Yield surfaces
#include "custom_constitutive/yield_surfaces/finite_strain/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/finite_strain/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/finite_strain/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/finite_strain/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/finite_strain/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/finite_strain/drucker_prager_plastic_potential.h"

namespace Kratos
{

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

/**
 * @class KratosStructuralMechanicsApplicationConstitutiveLaws
 * @ingroup StructuralMechanicsApplication
 * @brief This application features Elements, Conditions, Constitutive laws and Utilities for structural analysis problems
 * @author Alejandro Cornejo
 */
class  KratosStructuralMechanicsApplicationConstitutiveLaws : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosStructuralMechanicsApplicationConstitutiveLaws
    KRATOS_CLASS_POINTER_DEFINITION(KratosStructuralMechanicsApplicationConstitutiveLaws);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosStructuralMechanicsApplicationConstitutiveLaws();

    /// Destructor.
    ~KratosStructuralMechanicsApplicationConstitutiveLaws() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;

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
    std::string Info() const override
    {
        return "KratosStructuralMechanicsApplicationConstitutiveLaws";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }


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

    // Damage and plasticity laws
    const SmallStrainIsotropicPlasticityFactory mSmallStrainIsotropicPlasticityFactory;
    const FiniteStrainIsotropicPlasticityFactory mFiniteStrainIsotropicPlasticityFactory;
    const SmallStrainIsotropicDamageFactory mSmallStrainIsotropicDamageFactory;
    const ViscousGeneralizedKelvin<ElasticIsotropic3D> mViscousGeneralizedKelvin3D;
    const ViscousGeneralizedMaxwell<ElasticIsotropic3D> mViscousGeneralizedMaxwell3D;
    const GenericSmallStrainViscoplasticity3D mGenericSmallStrainViscoplasticity3D;

    /// Plasticity
    /* Small strain */
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerTresca;

    /* Finite strain */
    // Kirchhoff
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerTresca;

    // Neo-Hookean
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerTresca;

    /// Damage
    /* Small strain */
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuTresca;

    // d+d- laws
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D;
    
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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosStructuralMechanicsApplicationConstitutiveLaws& operator=(KratosStructuralMechanicsApplicationConstitutiveLaws const& rOther);

    /// Copy constructor.
    KratosStructuralMechanicsApplicationConstitutiveLaws(KratosStructuralMechanicsApplicationConstitutiveLaws const& rOther);


};// Class KratosStructuralMechanicsApplicationConstitutiveLaws

}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CONSTITUTIVE_LAWS_H_INCLUDED  defined