from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.HeatEquationApplication
#import KratosMultiphysics.MeshingApplication as KratosMeshing

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return HeatEquationSolver(main_model_part, custom_settings)

class HeatEquationSolver(object):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        base_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "heat_equation_solver",
            "problem_data" :{
                "domain_size" : 2
                },

            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
                },
            "relative_tolerance"        : 1e-5,
            "absolute_tolerance"        : 1e-7,
            "time_order"                : 2,
            "echo_level"                : 0,
            "maximum_iterations"        : 20,
            "compute_reactions"         : false,
            "reform_dofs_at_each_step"  : false,
            "move_mesh_flag"            : false,
            "domain_model_part"         : "your_computational_domain",
            "dynamic_tau"               : 1
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(base_settings)
        # self.CalculateReactionFlag = True
        ## Construct the linear solver

        #import linear_solver_factory
        #self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
        self.linear_solver = SuperLUIterativeSolver()

        print("Construction of HeatEquationSolver finished")

        # self.CalculateNodalAreaProcess = CalculateNodalAreaProcess(main_model_part, DOMAIN_SIZE)

    def AddVariables(self):
        ## Add base class variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY);
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY);
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT);
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX);
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE);
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY);
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)  ## check if it is ok
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX )
        print("Base class Heat Equation solver variable added correctly")

    def ImportModelPart(self):
        ## Read model part
        self._ModelPartReading()
        ## Replace elements and conditions
        self._ExecuteAfterReading()
        ## Set buffer size
        self._SetBufferSize()

        # print(self.main_model_part.GetSubModelPart("Wall_solid"))
        print ("Base class model reading finished.")

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()

        ## Model part writing
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def AddDofs(self):
        ## Adding dofs
        for node in self.main_model_part.Nodes:
            node.AddDof(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.FACE_HEAT_FLUX)

        print("Base class Heat Equation solver DOFs added correctly.")

    def AdaptMesh(self):
        pass

    def GetMinimumBufferSize(self):
        return 3

    def GetComputingModelPart(self):
        return self.main_model_part #.GetSubModelPart(self.settings["domain_model_part"].GetString())

    def GetOutputVariables(self):
        pass

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # Creating the solution strategy
        self.conv_criteria = KratosMultiphysics.ResidualCriteria(self.settings["relative_tolerance"].GetDouble(),
                                                                 self.settings["absolute_tolerance"].GetDouble())

        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,
                                                                            self.settings["time_order"].GetInt())

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        ## change 2 to DOMAIN_SIZE
        ## KratosMultiphysics.CalculateNodalAreaProcess(self.computing_model_part, DOMAIN_SIZE).Execute()
        # KratosMultiphysics.CalculateNodalAreaProcess(self.computing_model_part, 2).Execute()

        KratosMultiphysics.CalculateNodalAreaProcess(self.computing_model_part,
                                                     self.settings["problem_data"]["domain_size"].GetInt()).Execute()
        # for elem in self.main_model_part.Elements:
        #     print(elem)
        # err
        (self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.
        (self.solver).Check()


    def SaveRestart(self):
        pass #one should write the restart file here

    def Clear(self):
        (self.solver).Clear()

    def Check(self):
        (self.solver).Check()

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def SolverInitialize(self):
        (self.solver).Initialize()

    def SolverInitializeSolutionStep(self):
        (self.bdf_process).Execute()
        (self.solver).InitializeSolutionStep()

    def SolverPredict(self):
        (self.solver).Predict()

    def SolverSolveSolutionStep(self):
        (self.solver).SolveSolutionStep()

    def SolverFinalizeSolutionStep(self):
        (self.solver).FinalizeSolutionStep()

    def Solve(self):
        #self.SolverInitialize()

        # for elem in self.main_model_part.Elements:
            # print(elem)
        # err
        self.SolverInitializeSolutionStep()

        self.SolverPredict()

        self.SolverSolveSolutionStep()

        self.SolverFinalizeSolutionStep()

    def _ModelPartReading(self):
        ## Model part reading
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            ## Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
        else:
            raise Exception("Other input options are not yet implemented.")

    def _ExecuteAfterReading(self):
        # Read the DENSITY, CONDUCTIVITY and SPECIFIC_HEAT we apply it to the nodes
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            k = el.Properties.GetValue(KratosMultiphysics.CONDUCTIVITY)
            cp = el.Properties.GetValue(KratosMultiphysics.SPECIFIC_HEAT)
            break

        print(self.settings["dynamic_tau"].GetDouble())
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TAU] = self.settings["dynamic_tau"].GetDouble()

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.CONDUCTIVITY, k, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.SPECIFIC_HEAT, cp, self.main_model_part.Nodes)

    def _SetBufferSize(self):
        ## Set the buffer size
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize(self.GetMinimumBufferSize())
