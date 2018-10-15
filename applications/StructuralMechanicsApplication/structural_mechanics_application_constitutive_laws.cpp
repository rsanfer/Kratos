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

// System includes

// External includes

// Project includes
#include "structural_mechanics_application_constitutive_laws.h"

namespace Kratos {

    KratosStructuralMechanicsApplicationConstitutiveLaws::KratosStructuralMechanicsApplicationConstitutiveLaws()
        : KratosApplication("KratosStructuralMechanicsApplicationConstitutiveLaws")
        {}

    void KratosStructuralMechanicsApplicationConstitutiveLaws::Register() {
        // calling base class register to register Kratos components
        KratosApplication::Register();

        // Damage and plasticity
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticityFactory", mSmallStrainIsotropicPlasticityFactory);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamageFactory", mSmallStrainIsotropicDamageFactory);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("ViscousGeneralizedKelvin3D", mViscousGeneralizedKelvin3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("ViscousGeneralizedMaxwell3D", mViscousGeneralizedMaxwell3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("GenericSmallStrainViscoplasticity3D", mGenericSmallStrainViscoplasticity3D);

        // Custom Constitutive laws
        /// Plasticity

        /* Small strain */
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesVonMises", mSmallStrainIsotropicPlasticity3DVonMisesVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesDruckerPrager", mSmallStrainIsotropicPlasticity3DVonMisesDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesTresca", mSmallStrainIsotropicPlasticity3DVonMisesTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaVonMises", mSmallStrainIsotropicPlasticity3DTrescaVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaDruckerPrager", mSmallStrainIsotropicPlasticity3DTrescaDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaTresca", mSmallStrainIsotropicPlasticity3DTrescaTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerVonMises", mSmallStrainIsotropicPlasticity3DDruckerPragerVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager", mSmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerTresca", mSmallStrainIsotropicPlasticity3DDruckerPragerTresca);

        /* Finite strain */

        // Kirchhoff
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesModifiedMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesTresca", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombTresca", mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaModifiedMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaTresca", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerModifiedMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerTresca", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerTresca);

        // Neo-Hookean
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesModifiedMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaModifiedMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerModifiedMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerTresca);

        /// Damage
        /* Small strain */
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesVonMises", mSmallStrainIsotropicDamage3DVonMisesVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesDruckerPrager", mSmallStrainIsotropicDamage3DVonMisesDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesTresca", mSmallStrainIsotropicDamage3DVonMisesTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises", mSmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombTresca", mSmallStrainIsotropicDamage3DModifiedMohrCoulombTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaVonMises", mSmallStrainIsotropicDamage3DTrescaVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaDruckerPrager", mSmallStrainIsotropicDamage3DTrescaDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaTresca", mSmallStrainIsotropicDamage3DTrescaTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerVonMises", mSmallStrainIsotropicDamage3DDruckerPragerVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerDruckerPrager", mSmallStrainIsotropicDamage3DDruckerPragerDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerTresca", mSmallStrainIsotropicDamage3DDruckerPragerTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineVonMises", mSmallStrainIsotropicDamage3DRankineVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineDruckerPrager", mSmallStrainIsotropicDamage3DRankineDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineTresca", mSmallStrainIsotropicDamage3DRankineTresca);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuVonMises", mSmallStrainIsotropicDamage3DSimoJuVonMises);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuDruckerPrager", mSmallStrainIsotropicDamage3DSimoJuDruckerPrager);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuTresca", mSmallStrainIsotropicDamage3DSimoJuTresca);

        // d+d- laws
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineRankine3D", mSmallStrainDplusDminusDamageRankineRankine3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineSimoJu3D", mSmallStrainDplusDminusDamageRankineSimoJu3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineVonMises3D", mSmallStrainDplusDminusDamageRankineVonMises3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineTresca3D", mSmallStrainDplusDminusDamageRankineTresca3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineDruckerPrager3D", mSmallStrainDplusDminusDamageRankineDruckerPrager3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuRankine3D", mSmallStrainDplusDminusDamageSimoJuRankine3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuSimoJu3D", mSmallStrainDplusDminusDamageSimoJuSimoJu3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuVonMises3D", mSmallStrainDplusDminusDamageSimoJuVonMises3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuTresca3D", mSmallStrainDplusDminusDamageSimoJuTresca3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuDruckerPrager3D", mSmallStrainDplusDminusDamageSimoJuDruckerPrager3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesRankine3D", mSmallStrainDplusDminusDamageVonMisesRankine3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesSimoJu3D", mSmallStrainDplusDminusDamageVonMisesSimoJu3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesVonMises3D", mSmallStrainDplusDminusDamageVonMisesVonMises3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesTresca3D", mSmallStrainDplusDminusDamageVonMisesTresca3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesDruckerPrager3D", mSmallStrainDplusDminusDamageVonMisesDruckerPrager3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaRankine3D", mSmallStrainDplusDminusDamageTrescaRankine3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaSimoJu3D", mSmallStrainDplusDminusDamageTrescaSimoJu3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaVonMises3D", mSmallStrainDplusDminusDamageTrescaVonMises3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaTresca3D", mSmallStrainDplusDminusDamageTrescaTresca3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaDruckerPrager3D", mSmallStrainDplusDminusDamageTrescaDruckerPrager3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerRankine3D", mSmallStrainDplusDminusDamageDruckerPragerRankine3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerSimoJu3D", mSmallStrainDplusDminusDamageDruckerPragerSimoJu3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerVonMises3D", mSmallStrainDplusDminusDamageDruckerPragerVonMises3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerTresca3D", mSmallStrainDplusDminusDamageDruckerPragerTresca3D);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D", mSmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D);

    }
} // namespace Kratos.