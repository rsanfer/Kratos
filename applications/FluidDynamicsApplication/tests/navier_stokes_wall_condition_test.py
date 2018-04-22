import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

import KratosMultiphysics.KratosUnittest as UnitTest

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class NavierStokesWallConditionTest(UnitTest.TestCase):
    def testNavierStokesWallCondition(self):
        self.setUp()
        self.setUpProblem()
        self.runTest()
        self.tearDown()
        self.checkResults()

    def setUp(self):
        self.ext_pres_res = 1000.0
        self.check_tolerance = 2e-2
        self.print_output = False
        self.work_folder = "NavierStokesWallConditionTest"   
        self.settings = "NavierStokesWallConditionTestParameters.json"

    def tearDown(self):
        with WorkFolderScope(self.work_folder):
            try:
                os.remove(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')
            except FileNotFoundError as e:
                pass

    def setUpProblem(self):
        with WorkFolderScope(self.work_folder):
            with open(self.settings, 'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

            Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

            ## Solver construction
            import python_solvers_wrapper_fluid
            self.solver = python_solvers_wrapper_fluid.CreateSolver(self.main_model_part, self.ProjectParameters)

            self.solver.AddVariables()

            ## Read the model - note that SetBufferSize is done here
            self.solver.ImportModelPart()

            ## Add AddDofs
            self.solver.AddDofs()

            ## Solver initialization
            self.solver.Initialize()

            ## Get the list of the skin submodel parts in the object Model
            for i in range(self.ProjectParameters["solver_settings"]["skin_parts"].size()):
                skin_part_name = self.ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
                Model.update({skin_part_name: self.main_model_part.GetSubModelPart(skin_part_name)})

            ## Get the gravity submodel part in the object Model
            for i in range(self.ProjectParameters["gravity"].size()):
                gravity_part_name = self.ProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
                Model.update({gravity_part_name: self.main_model_part.GetSubModelPart(gravity_part_name)})

            ## Processes construction
            import process_factory
            self.list_of_processes  = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( self.ProjectParameters["gravity"] )
            self.list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( self.ProjectParameters["boundary_conditions_process_list"] )

            ## Processes initialization
            for process in self.list_of_processes:
                process.ExecuteInitialize()

    def runTest(self):
        with WorkFolderScope(self.work_folder):
            if (self.print_output):
                gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
                multifile = KratosMultiphysics.MultiFileFlag.SingleFile
                deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
                write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
                gid_io = KratosMultiphysics.GidIO(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString(),gid_mode,multifile,deformed_mesh_flag, write_conditions)

                mesh_name = 0.0
                gid_io.InitializeMesh( mesh_name)
                gid_io.WriteMesh( self.main_model_part.GetMesh() )
                gid_io.FinalizeMesh()
                gid_io.InitializeResults(mesh_name,(self.main_model_part).GetMesh())

            end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

            time = 0.0
            step = 0

            for process in self.list_of_processes:
                process.ExecuteBeforeSolutionLoop()

            while(time <= end_time):

                Dt = self.solver.ComputeDeltaTime()
                step += 1
                time += Dt
                self.main_model_part.CloneTimeStep(time)
                self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = step

                for process in self.list_of_processes:
                    process.ExecuteInitializeSolutionStep()

                self.solver.Solve()

                for process in self.list_of_processes:
                    process.ExecuteFinalizeSolutionStep()

                for process in self.list_of_processes:
                    process.ExecuteBeforeOutputStep()

                if (self.print_output):
                    gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY,self.main_model_part.Nodes,time,0)
                    gid_io.WriteNodalResults(KratosMultiphysics.PRESSURE,self.main_model_part.Nodes,time,0)
                    gid_io.WriteNodalResults(KratosMultiphysics.DISTANCE,self.main_model_part.Nodes,time,0)

                for process in self.list_of_processes:
                    process.ExecuteAfterOutputStep()

            for process in self.list_of_processes:
                process.ExecuteFinalize()

            if (self.print_output):
                gid_io.FinalizeResults()

    def checkResults(self):
        react_x = 0.0

        for node in self.main_model_part.Nodes:
            react_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)

        self.assertAlmostEqual(-react_x, self.ext_pres_res, delta=2e-2)

if __name__ == '__main__':
    test = NavierStokesWallConditionTest()
    test.setUp()
    test.setUpProblem()
    test.runTest()
    test.tearDown()
    test.checkResults()
    