
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.HeatEquationApplication as KratosHEq
KratosMultiphysics.CheckForPreviousImport()

import navier_stokes_base_solver
def CreateSolver(fluid_main_model_part, heat_equation_main_model_part, project_parameters):
         return CFDThermalSolver(fluid_main_model_part, heat_equation_main_model_part, project_parameters)

class CFDThermalSolver:
    """          """
    def __init__(self, fluid_main_model_part, heat_equation_main_model_part, project_parameters ):

        # Initial tests

        #obtain the compute_model_part from the MODEL once the object is implemented

        self.fluid_main_model_part = fluid_main_model_part
        self.heat_equation_main_model_part = heat_equation_main_model_part
        # print(self.fluid_main_model_part)
        # print(self.heat_equation_main_model_part)

        # Settings string in JSON format
        solvers_default_settings = KratosMultiphysics.Parameters("""
        {
        "heat_equation_solver_settings":
            {
            "solver_type": "heat_equation_solver",
            "problem_data" :{
                "domain_size" : 2
            },
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "coupled_solid"
            },
            "relative_tolerance"        : 1e-5,
            "absolute_tolerance"        : 1e-7,
            "time_order"                : 2,
            "echo_level"                : 0,
            "maximum_iterations"        : 20,
            "compute_reactions"         : true,
            "reform_dofs_at_each_step"  : false,
            "move_mesh_flag"            : false,
            "domain_model_part"         : "your_computational_domain",
            "dynamic_tau"               : 1
            },
        "fluid_solver_settings":
            {
            "solver_type": "navier_stokes_solver_vmsmonolithic",
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "coupled_fluid_vel"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.0,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "compute_reactions": true,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_step": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        :
            {
                    "solver_type" : "AMGCL",
                    "smoother_type":"ilu0",
                    "krylov_type": "gmres",
                    "coarsening_type": "aggregation",
                    "max_iteration": 200,
                    "provide_coordinates": false,
                    "gmres_krylov_space_dimension": 100,
                    "verbosity" : 0,
                    "tolerance": 1e-7,
                    "scaling": false,
                    "block_size": 1,
                    "use_block_matrices_if_possible" : true,
                    "coarse_enough" : 5000
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                :
            {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "user_delta_time"     : 0.004
            },
            "alpha":-0.3,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "move_mesh_flag": false,
            "turbulence_model": "None"
            }
        }""")


        # Take the each one of the solvers settings from the ProjectParameters
        self.settings = KratosMultiphysics.Parameters("{}")
        self.settings.AddValue("heat_equation_solver_settings",project_parameters["heat_settings"]["heat_equation_solver_settings"])
        self.settings.AddValue("fluid_solver_settings",project_parameters["fluid_settings"]["fluid_solver_settings"])

        # Overwrite the default settings with user-provided parameters
        self.settings.RecursivelyValidateAndAssignDefaults(solvers_default_settings)
        # Auxiliar variables

        print("*** Partitioned CFDThermalSolver  construction starts...")
        # Construct the fluid solver
        fluid_solver_module = __import__(self.settings["fluid_solver_settings"]["solver_type"].GetString())
        self.fluid_solver = fluid_solver_module.CreateSolver(self.fluid_main_model_part,
                                                             self.settings["fluid_solver_settings"])

        print("* Fluid solver constructed.")

        # Construct the heat_equation solver
        heat_equation_solver_module = __import__(self.settings["heat_equation_solver_settings"]["solver_type"].GetString())
        self.heat_equation_solver = heat_equation_solver_module.CreateSolver(self.heat_equation_main_model_part,
                                                                             self.settings["heat_equation_solver_settings"])

        print("* Heat solver constructed.")


    def GetMinimumBufferSize(self):
        # Get fluid buffer size
        buffer_fluid = self.fluid_solver.GetMinimumBufferSize()
        # Get heat equation buffer size
        buffer_heat = self.heat_equation_solver.GetMinimumBufferSize()

        return min(buffer_heat,buffer_fluid)

    def AddVariables(self):

        self.fluid_solver.AddVariables()
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.heat_equation_solver.AddVariables()


    def ImportModelPart(self):
        #read fluid
        self.fluid_solver.ImportModelPart()
        # print(self.fluid_main_model_part)

        #read solid
        # self.heat_equation_so lver.ImportModelPart()
        # print(self.heat_equation_main_model_part)
        #print(self.heat_equation_solver.main_model_part.GetSubModelPart("Wall_solid"))
        # print(self.fluid_solver.GetComputingModelPart())


        #construct thermal model part from fluid (ConnectivityPreserveModelPart)
        #self.heat_equation_solver.ImportModelPart()
        #print(self.heat_equation_solver.GetComputingModelPart())

        #print(dir(self.fluid_solver))
        #print(help(self.fluid_solver.GetComputingModelPart()))
        mpfluid = self.fluid_solver.main_model_part #GetComputingModelPart()
        # print(mpfluid)

        mpheat = self.heat_equation_solver.main_model_part
        # print(mpheat)





        Dim =self.settings["heat_equation_solver_settings"]["problem_data"]["domain_size"].GetInt()
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if (Dim == 2):
            modeler.GenerateModelPart(mpfluid,
                                      mpheat,
                                      "HeatEquation2D", "LineCondition2D2N")
        else:
            modeler.GenerateModelPart(mpfluid,
                                      mpheat,
                                     "HeatEquation3D", "WallCondition3D3N")



        #print(self.heat_equation_solver.main_model_part.GetSubModelPart("Wall_solid"))
        print("********ConnectivityPreserveModeler done!***********")

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, 1.225, self.heat_equation_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.CONDUCTIVITY, 0.0257, self.heat_equation_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.SPECIFIC_HEAT, 1.005, self.heat_equation_solver.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.HEAT_FLUX, 0, self.heat_equation_solver.main_model_part.Nodes)



    def AddDofs(self):

        # Add DOFs fluid
        self.fluid_solver.AddDofs()
        # Add DOFs heat_equation
        self.heat_equation_solver.AddDofs()


    def Initialize(self):
        # Initialize fluid solver
        self.fluid_solver.Initialize()
        # Initialize heat_equation solver
        self.heat_equation_solver.Initialize()
    # #
    def SolverInitializeSolutionStep(self):
        # (self.bdf_process).Execute()
        self.fluid_solver.InitializeSolutionStep()
        self.heat_equation_solver.SolverInitializeSolutionStep()
    #
    def SolverPredict(self):
        self.fluid_solver.Predict()
        self.heat_equation_solver.SolverPredict()

    def SolverSolveSolutionStep(self):
        self.fluid_solver.SolveSolutionStep()
        self.heat_equation_solver.SolverSolveSolutionStep()

    def SolverFinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.heat_equation_solver.SolverFinalizeSolutionStep()

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()


    def GetOutputVariables(self):
        pass



    def SaveRestart(self):
        pass


    def Solve(self):
        # self.buoyancy_process.ExecuteInitializeSolutionStep()
        print("------------------fluid solving")
        self.fluid_solver.Solve()
        print("------------------Thermal solving")
        # for elem in self.heat_equation_solver.main_model_part.Elements:
        #     print(elem)
        # err
        self.heat_equation_solver.Solve()
