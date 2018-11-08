//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
//

#if !defined(KRATOS_FEM_TO_DEM_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_FEM_TO_DEM_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
	KRATOS_DEFINE_VARIABLE( double, DAMAGE_EDGE1)
	KRATOS_DEFINE_VARIABLE( double, DAMAGE_EDGE2)
	KRATOS_DEFINE_VARIABLE( double, DAMAGE_EDGE3)
	KRATOS_DEFINE_VARIABLE( double, DAMAGE_ELEMENT)
	KRATOS_DEFINE_VARIABLE(Vector, STRESS_VECTOR);
	KRATOS_DEFINE_VARIABLE(double, YIELD_STRESS_C);
	KRATOS_DEFINE_VARIABLE(double, YIELD_STRESS_T);
	KRATOS_DEFINE_VARIABLE(int, ITER);
	KRATOS_DEFINE_VARIABLE(double, FRAC_ENERGY_T)
	KRATOS_DEFINE_VARIABLE(double, FRAC_ENERGY_C)
	KRATOS_DEFINE_VARIABLE(double, NODAL_DAMAGE)
	KRATOS_DEFINE_VARIABLE(Vector, STRESS_VECTOR_INTEGRATED);
	KRATOS_DEFINE_VARIABLE(double, THRESHOLD)
	KRATOS_DEFINE_VARIABLE(Vector, SMOOTHED_STRESS_VECTOR);
	KRATOS_DEFINE_VARIABLE(std::string, YIELD_SURFACE);
	KRATOS_DEFINE_VARIABLE(Vector, STRAIN_VECTOR);
	KRATOS_DEFINE_VARIABLE(bool, TANGENT_CONSTITUTIVE_TENSOR);
	KRATOS_DEFINE_VARIABLE(bool, SMOOTHING);
	KRATOS_DEFINE_VARIABLE(double, IS_DAMAGED);
	KRATOS_DEFINE_VARIABLE(int, MESH_REFINED);
	KRATOS_DEFINE_VARIABLE(int, IS_DYNAMIC);
	KRATOS_DEFINE_VARIABLE(double, STRESS_THRESHOLD);
	KRATOS_DEFINE_VARIABLE(double, INITIAL_THRESHOLD);
	KRATOS_DEFINE_VARIABLE(int, INTEGRATION_COEFFICIENT);
	KRATOS_DEFINE_VARIABLE(std::string, MAPPING_PROCEDURE);
	KRATOS_DEFINE_VARIABLE(bool, IS_DEM);
	KRATOS_DEFINE_VARIABLE(double, DEM_RADIUS);
	KRATOS_DEFINE_VARIABLE(bool, DEM_GENERATED);
	KRATOS_DEFINE_VARIABLE(bool, INACTIVE_NODE);
	KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_ACTIVE_ELEMENTS);
	KRATOS_DEFINE_VARIABLE(bool, NODAL_FORCE_APPLIED);
	KRATOS_DEFINE_VARIABLE(double, NODAL_FORCE_X);
	KRATOS_DEFINE_VARIABLE(double, NODAL_FORCE_Y);
	KRATOS_DEFINE_VARIABLE(double, NODAL_FORCE_Z);
	KRATOS_DEFINE_VARIABLE(Vector, NODAL_STRESS_VECTOR);
	KRATOS_DEFINE_VARIABLE(double, EQUIVALENT_NODAL_STRESS);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(EQUIVALENT_NODAL_STRESS_GRADIENT);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUXILIAR_GRADIENT);
	

	KRATOS_DEFINE_VARIABLE(Matrix, STRAIN_TENSOR);
	KRATOS_DEFINE_VARIABLE(Matrix, STRESS_TENSOR);
	KRATOS_DEFINE_VARIABLE(Matrix, STRESS_TENSOR_INTEGRATED);
	
	// composite
	KRATOS_DEFINE_VARIABLE(Matrix, CONCRETE_STRESS_TENSOR);
	KRATOS_DEFINE_VARIABLE(Matrix, STEEL_STRESS_TENSOR);
	KRATOS_DEFINE_VARIABLE(Vector, CONCRETE_STRESS_VECTOR);
	KRATOS_DEFINE_VARIABLE(Vector, STEEL_STRESS_VECTOR);
	KRATOS_DEFINE_VARIABLE(double,YOUNG_MODULUS_STEEL);
	KRATOS_DEFINE_VARIABLE(double,DENSITY_STEEL);
	KRATOS_DEFINE_VARIABLE(double,POISSON_RATIO_STEEL);
	KRATOS_DEFINE_VARIABLE(double,STEEL_VOLUMETRIC_PART);
	KRATOS_DEFINE_VARIABLE(Matrix,CONCRETE_STRESS_TENSOR_INTEGRATED);
	
	//plasticity steel
	KRATOS_DEFINE_VARIABLE(double,YIELD_STRESS_C_STEEL);
	KRATOS_DEFINE_VARIABLE(double,YIELD_STRESS_T_STEEL);
	KRATOS_DEFINE_VARIABLE(double,FRACTURE_ENERGY_STEEL);
	KRATOS_DEFINE_VARIABLE(double,PLASTIC_DISSIPATION_CAPAP);
	KRATOS_DEFINE_VARIABLE(double,EQUIVALENT_STRESS_VM);
	KRATOS_DEFINE_VARIABLE(int,HARDENING_LAW);
	KRATOS_DEFINE_VARIABLE(double,MAXIMUM_STRESS);
	KRATOS_DEFINE_VARIABLE(double,MAXIMUM_STRESS_POSITION);
	KRATOS_DEFINE_VARIABLE(bool, IS_TAKEN);
	KRATOS_DEFINE_VARIABLE(int, PRESSURE_ID);
	
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD);
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD);
	
	
	
}

#endif	/* KRATOS_FEM_TO_DEM_APPLICATION_VARIABLES_H_INCLUDED */
