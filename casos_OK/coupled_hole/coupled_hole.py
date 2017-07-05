##########################################################
##################################################################
#setting the domain size for the problem to be solved

#including kratos path
import sys
from sympy import *
import math
import numpy as np
import matplotlib.pyplot as plt
from KratosMultiphysics import *    #we import the KRATOS
from KratosMultiphysics.HeatEquationApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem
import KratosMultiphysics.ExternalSolversApplication



#defining a model part
# model_part = ModelPart("ExampleModelPart");  #we create a model part

ProjectParameters = Parameters("""
{
"heat_settings" :
{
"problem_data" :{
    "domain_size" : 2,
    "model_part_name" :"thermal_main_model_part"
},
"heat_equation_solver_settings": {
    "solver_type": "heat_equation_solver",
    "model_import_settings": {
        "input_type": "mdpa",
        "input_filename": "coupled_fluid_temp"
    },
    "relative_tolerance"        : 1e-5,
    "absolute_tolerance"        : 1e-7,
    "time_order"                : 2,
    "echo_level"                : 0,
    "maximum_iterations"        : 20,
    "compute_reactions"         : false,
    "reform_dofs_at_each_step"  : false,
    "move_mesh_flag"            : false,
    "domain_model_part"         : "Parts_Parts_Auto1",
    "dynamic_tau"               : 1
}
},
"fluid_settings" :
{
"problem_data" : {
    "domain_size" : 2,
    "model_part_name" :"fluid_main_model_part"
},
"fluid_solver_settings": {
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
    "compute_reactions": false,
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
    "volume_model_part_name" : "Parts_Parts_Auto1",
    "skin_parts": ["AutomaticInlet2D_Automatic_inlet_velocity_Auto2","NoSlip2D_walls","NoSlip2D_No_Slip_Cylinder", "VelocityConstraints2D_outlet"],
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
}
}""")
## thermo-fluid model parts definition
fluid_main_model_part = ModelPart(ProjectParameters["fluid_settings"]["problem_data"]["model_part_name"].GetString())
fluid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["fluid_settings"]["problem_data"]["domain_size"].GetInt())

thermal_main_model_part = ModelPart(ProjectParameters["heat_settings"]["problem_data"]["model_part_name"].GetString())
thermal_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["heat_settings"]["problem_data"]["domain_size"].GetInt())




# FluidModel = {ProjectParameters["fluid_settings"]["problem_data"]["model_part_name"].GetString() : fluid_main_model_part}
# ThermalModel = {ProjectParameters["heat_settings"]["problem_data"]["model_part_name"].GetString() : thermal_main_model_part}



# import CFD_Thermal_solver           #we import the python file that includes the commands that we need
# heat_equation_solver = heat_equation_solver.CreateSolver(model_part, settings)

solver_module = __import__("CFD_thermal_solver")
solver = solver_module.CreateSolver(fluid_main_model_part,thermal_main_model_part, ProjectParameters)




# thermal_main_model_part.AddNodalSolutionStepVariable(BODY_FORCE)
# for node in thermal_main_model_part.Nodes:
#     node.AddDof(BODY_FORCE)
# fluid_main_model_part.AddNodalSolutionStepVariable(BODY_FORCE)
# for node in fluid_main_model_part.Nodes:
#     node.AddDof(BODY_FORCE)

solver.AddVariables()

solver.ImportModelPart()



#now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file
gid_mode = GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("art4",gid_mode,multifile,deformed_mesh_flag,write_conditions)


# we create a mesh for the postprocess

out_mesh = fluid_main_model_part.GetMesh()

out_nodes = fluid_main_model_part.Nodes

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh(out_mesh);

gid_io.FinalizeMesh()


# we add the DoFs
solver.AddDofs()


solver.Initialize()


#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)
fluid_main_model_part.SetBufferSize(3)
thermal_main_model_part.SetBufferSize(3)

# Initial conditions
# Temperature fluid initial condition
initial_fluid_temperature = 298.15
#Temperature solid initial condition
initial_solid_temperature = 355.15


time = 0.0
delta_time = 0.05

total_steps = 100
v_avarage = 1.5
# CalculateNodalAreaProcess(fluid_main_model_part,ProjectParameters["heat_settings"]["problem_data"]["domain_size"].GetInt()).Execute()

#### fluid
for node in fluid_main_model_part.GetSubModelPart("AutomaticInlet2D_Automatic_inlet_velocity_Auto2").Nodes:
    Inlet_velocity = Vector(3)
    # inlet BC  U (0, y) = 4 * U_m * y(H − y)/H^2 , V = 0
    # H = 0.41
    Inlet_velocity[0] = 4*v_avarage*node.Y*(0.41-node.Y)/(0.41*0.41)
    Inlet_velocity[1] = 0.0
    Inlet_velocity[2] = 0.0
    node.SetSolutionStepValue(VELOCITY, Inlet_velocity)
    node.Fix(VELOCITY_X)
    node.Fix(VELOCITY_Y)
    node.Fix(VELOCITY_Z)
    #fluid Temperature at the inlet
for node in thermal_main_model_part.GetSubModelPart("AutomaticInlet2D_Automatic_inlet_velocity_Auto2").Nodes:
    node.SetSolutionStepValue(TEMPERATURE,initial_fluid_temperature)
    node.Fix(TEMPERATURE)

# outlet   V = 0
for node in fluid_main_model_part.GetSubModelPart("VelocityConstraints2D_outlet").Nodes:
    Outlet_velocity = Vector(3)
    Outlet_velocity[1] = 0.0
    node.SetSolutionStepValue(VELOCITY, Outlet_velocity)
    node.SetSolutionStepValue(PRESSURE, 0.0)
    node.Fix(VELOCITY_Y)
    node.Fix(PRESSURE)


fluid_wall_velocity = Vector(3)
fluid_wall_velocity[0] = 0.0
fluid_wall_velocity[1] = 0.0
fluid_wall_velocity[2] = 0.0
# no-slip condition
for node in fluid_main_model_part.GetSubModelPart("NoSlip2D_walls").Nodes:
    node.SetSolutionStepValue(VELOCITY, fluid_wall_velocity)
    node.Fix(VELOCITY_X)
    node.Fix(VELOCITY_Y)
    node.Fix(VELOCITY_Z)


for node in fluid_main_model_part.GetSubModelPart("NoSlip2D_No_Slip_Cylinder").Nodes:
    node.SetSolutionStepValue(VELOCITY, fluid_wall_velocity)
    node.Fix(VELOCITY_X)
    node.Fix(VELOCITY_Y)
    node.Fix(VELOCITY_Z)

for node in thermal_main_model_part.Nodes:
    node.SetSolutionStepValue(TEMPERATURE,initial_fluid_temperature)


for node in  thermal_main_model_part.GetSubModelPart("NoSlip2D_No_Slip_Cylinder").Nodes:
    node.SetSolutionStepValue(TEMPERATURE,initial_solid_temperature)
    # node.Fix(TEMPERATURE)

from time import sleep
rho0 = 1.22500
g0 = 9.81
for step in range(1,total_steps):

    time += delta_time
    fluid_main_model_part.CloneTimeStep(time)

    if (step >= 3):
        for node in thermal_main_model_part.Nodes:


            Temperature = node.GetSolutionStepValue(TEMPERATURE)
            rho = rho0 * (1 - 1/initial_fluid_temperature * (Temperature - initial_fluid_temperature))
            # g = g0 * (1 - 1/initial_fluid_temperature * (Temperature - initial_fluid_temperature))
            # print(rho)
            # print(g)
            node.SetSolutionStepValue(DENSITY, rho)
            # node.SetSolutionStepValue(BODY_FORCE, g)
        sleep(1)
        solver.Solve()
        print ("Solved!")

    gid_io.InitializeResults(time,out_mesh)
    # gid_io.WriteNodalResults(HEAT_FLUX,out_nodes,time,0)
    # gid_io.WriteNodalResults(NODAL_AREA,out_nodes,time,0)
    gid_io.WriteNodalResults(TEMPERATURE,out_nodes,time,0)
    gid_io.WriteNodalResults(VELOCITY,out_nodes,time,0)
    gid_io.WriteNodalResults(PRESSURE,out_nodes,time,0)
    gid_io.WriteNodalResults(NODAL_AREA,out_nodes,time,0)
    gid_io.WriteNodalResults(DENSITY,out_nodes,time,0)
    gid_io.FinalizeResults()
