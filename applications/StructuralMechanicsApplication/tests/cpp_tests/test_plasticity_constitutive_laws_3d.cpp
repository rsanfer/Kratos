// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"

// Application includes

// Integrator
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
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
// Constitutive law
#include "custom_constitutive/generic_small_strain_isotropic_plasticity_3d.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"


namespace Kratos
{
namespace Testing
{
    // We test the associated plasticity Constitutive laws...
        typedef Node < 3 > NodeType;
        
        typedef GenericSmallStrainIsotropicPlasticity3D
			<GenericConstitutiveLawIntegratorPlasticity
				<ModifiedMohrCoulombYieldSurface
					<ModifiedMohrCoulombPlasticPotential>>> MC;

		typedef GenericSmallStrainIsotropicPlasticity3D
			<GenericConstitutiveLawIntegratorPlasticity
				<VonMisesYieldSurface
					<VonMisesPlasticPotential>>> VM;

		typedef GenericSmallStrainIsotropicPlasticity3D
			<GenericConstitutiveLawIntegratorPlasticity
				<DruckerPragerYieldSurface
					<DruckerPragerPlasticPotential>>> DP;

		typedef GenericSmallStrainIsotropicPlasticity3D
			<GenericConstitutiveLawIntegratorPlasticity
				<TrescaYieldSurface
					<TrescaPlasticPotential>>> T;

    /** 
    * Check the correct calculation of the integrated stress with the CL's
    */
    KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressPlasticity, KratosStructuralMechanicsFastSuite)
    {
        ConstitutiveLaw::Parameters rValues;
		Properties rMaterialProperties;
        Vector rStressVector, rStrainVector;

        ModelPart TestMdpa("Main");

        NodeType::Pointer Node1 = TestMdpa.CreateNewNode(1, 0.0, 0.0, 0.0);
        NodeType::Pointer Node2 = TestMdpa.CreateNewNode(2, 1.0, 0.0, 0.0);
        NodeType::Pointer Node3 = TestMdpa.CreateNewNode(3, 0.0, 1.0, 0.0);
        NodeType::Pointer Node4 = TestMdpa.CreateNewNode(4, 0.0, 0.0, 1.0);

        Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(Node1,Node2,Node3,Node4);

        rStressVector = ZeroVector(6);
        rStrainVector = ZeroVector(6);
        rStrainVector[0] = 0.0;
        rStrainVector[1] = 0.0;
        rStrainVector[2] = -8.0e-5;
        rStrainVector[3] = 0.0;
        rStrainVector[4] = 0.0;
        rStrainVector[5] = -1.6941e-21;

        rMaterialProperties.SetValue(YOUNG_MODULUS, 210e9);
        rMaterialProperties.SetValue(POISSON_RATIO, 0.22);
        rMaterialProperties.SetValue(YIELD_STRESS_COMPRESSION, 3.0e6);
        rMaterialProperties.SetValue(YIELD_STRESS_TENSION, 3.0e6);
        rMaterialProperties.SetValue(FRICTION_ANGLE, 32.0);
        rMaterialProperties.SetValue(DILATANCY_ANGLE, 16.0);
        rMaterialProperties.SetValue(SOFTENING_TYPE, 1);
        rMaterialProperties.SetValue(FRACTURE_ENERGY, 1000.0);
        rMaterialProperties.SetValue(HARDENING_CURVE, 0);

        rValues.SetElementGeometry(Geom);
        rValues.SetMaterialProperties(rMaterialProperties);
        rValues.SetStrainVector(rStrainVector);
        rValues.SetStressVector(rStressVector);
        
		// Create the CL's
		MC MohrCoulombCL = MC();
		VM VonMisesCL = VM();
		DP DruckerPragerCL = DP();
		T TrescaCL = T();

        std::vector<double> MCres, VMres, DPres, Tres;
        MCres = { -9.07262e+06, -9.07262e+06, -1.18548e+07, 0.0, 0.0, -2.94576e-11 };
        VMres = { -9.09508e+06, -9.09508e+06, -1.18098e+07, 0.0, 0.0, -2.87441e-11 };
        DPres = { -5.40984e+06, -5.40984e+06, -1.91803e+07, 0.0, 0.0, -1.45804e-10 };
        Tres  = { -9.09508e+06, -9.09508e+06, -1.18098e+07, 0.0, 0.0, -2.87441e-11 };

        Vector TestMC, TestVM, TestDP, TestT;
        MohrCoulombCL.CalculateMaterialResponseCauchy(rValues);
        TestMC = rValues.GetStressVector();

		VonMisesCL.CalculateMaterialResponseCauchy(rValues);
        TestVM = rValues.GetStressVector();

		DruckerPragerCL.CalculateMaterialResponseCauchy(rValues);
        TestDP = rValues.GetStressVector();

		TrescaCL.CalculateMaterialResponseCauchy(rValues);
        TestT = rValues.GetStressVector();

        //Check the results
        for (int comp = 0; comp < 6; comp++) {
            KRATOS_CHECK_NEAR(MCres[comp], TestMC[comp], 1.0e-3);
            KRATOS_CHECK_NEAR(VMres[comp], TestVM[comp], 1.0e-3);
            KRATOS_CHECK_NEAR(DPres[comp], TestDP[comp], 1.0e-3);
            KRATOS_CHECK_NEAR(Tres[comp],  TestT[comp],  1.0e-3);
        }
    }
}
}
