import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.HeatEquationApplication as KratosHeat
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyBoussinesqProcess(Model, settings["Parameters"])

class ApplyBoussinesqProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                                   :"model_part",
            "initial_fluid_temperature"                         :298.15,
            "gravity_x"                                         :0.0,
            "gravity_y"                                         :9.81,
            "gravity_z"                                         :0.0,
            "fluid_thermal_expansion"                           :""
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.model_part_name = Model[settings["model_part_name"].GetString()]
        self.ComputeBodyForceProcess = KratosFluid.ComputeBodyForceProcess(self.model_part_name, settings)


    def ExecuteInitializeSolutionStep(self):
        self.ComputeBodyForceProcess.ExecuteInitializeSolutionStep()

# for node in thermal_main_model_part.Nodes:
#
#     Temperature = node.GetSolutionStepValue(TEMPERATURE)
#     # rho = rho0 * (1 - 1/initial_fluid_temperature * (Temperature - initial_fluid_temperature))
#     g = g0 * (1 - 1/initial_fluid_temperature * (Temperature - initial_fluid_temperature))
#     # print(rho)
#     print(g)
#     gravity[1] = g
#     node.SetSolutionStepValue(BODY_FORCE, gravity)
#     node.Fix(BODY_FORCE_X)
#     node.Fix(BODY_FORCE_Y)
#     node.Fix(BODY_FORCE_Z)
