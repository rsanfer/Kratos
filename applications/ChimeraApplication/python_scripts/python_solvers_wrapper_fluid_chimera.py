from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_type = custom_settings["solver_settings"]["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "monolithic" or solver_type == "Monolithic"):
            solver_module_name = "navier_stokes_solver_vmsmonolithic_chimera"

        elif (solver_type == "fractional_step" or solver_type == "FractionalStep"):
            solver_module_name = "navier_stokes_solver_fractionalstep_chimera"
        else:
            raise Exception("the requested solver type is not in the python solvers wrapper")

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        raise Exception("the requested solver type is not in the python solvers wrapper")
    else:
        raise Exception("parallelism is neither OpenMP nor MPI")

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(model, custom_settings["solver_settings"])

    return solver
