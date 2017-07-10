     ##########################################################
##################################################################
#setting the domain size for the problem to be solved

#including kratos path
import sys
from time import sleep
from sympy import *
import math
import numpy as np
from KratosMultiphysics import *    #we import the KRATOS
from KratosMultiphysics.HeatEquationApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.FSIApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem
from KratosMultiphysics.FSIApplication import *
import NonConformant_OneSideMap

def ConnectivityMapper(fluid_nodes, solid_nodes):
    if (len(fluid_nodes)!=len(solid_nodes)):
        sys.exit("ERROR! The solid and the fluid have a mismatching nodes mesh at the interface!")
    else:
        print("Solid-Fluid interface matching nodes")

    connectivity_list=[]
    id_list=[]
    tol=1e-8 #Geometric tolerance.
    i=0

    for node in fluid_nodes:
        x=node.X
        y=node.Y
        suptolx=x+tol
        inftolx=x-tol
        suptoly=y+tol
        inftoly=y-tol
        #~ print('Para nodo con coordenadas ',x,' ; ',y)
        #~ print(suptolx,' ',inftolx,' ',suptoly,' ',inftoly)
        j=0
        for node_solid in solid_nodes:
            #~ print(node_solid.Id)
            #~ err
            x_solid=node_solid.X
            y_solid=node_solid.Y

            if x_solid<suptolx and x_solid>inftolx and y_solid<suptoly and y_solid>inftoly:
                #~ print('Tenemos un caso!')
                connectivity_list.append([i,j])             # List w/ positions: i position in fluid corresponds to j position in solid
                id_list.append([node.Id,node_solid.Id])     # List w/ nodal id's: same concept but with nodal Id's.
                print('Fluid node: ',node.Id,'; Solid node: ',node_solid.Id)
                #~ time.sleep(3)
                break

            j+=1

        i+=1
    return id_list
    print('Connectivity Mapper interface fluid_nodes/solid_nodes done!')





#defining a model e
# model_part = ModelPart("ExampleModelPart");  #we create a model part
ProjectParameters1 = Parameters("""
{
    "solver_type": "heat_equation_solver",
    "problem_data" :{
        "domain_size" : 2
        },
    "model_import_settings" : {
        "input_type": "mdpa",
        "input_filename": "cylinders"
        },
    "time_order" : 2,
    "echo_level" : 0,
    "domain_model_part": "Parts_COPPER"
}""")

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
                "input_filename": "HEturb_vel"
            },
            "relative_tolerance"        : 1e-5,
            "absolute_tolerance"        : 1e-7,
            "time_order"                : 2,
            "echo_level"                : 0,
            "maximum_iterations"        : 20,
            "compute_reactions"         : true,
            "reform_dofs_at_each_step"  : false,
            "move_mesh_flag"            : false,
            "domain_model_part"         : "Parts_AIR",
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
                "input_filename": "HEturb_vel"
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
            "volume_model_part_name" : "Parts_AIR",
            "skin_parts": ["AutomaticInlet2D_INLET","NoSlip2D_WALLS", "Outlet2D_OUTLET", "NoSlip2D_FLUID_INTERFACE"],
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

solid_main_model_part =  ModelPart(ProjectParameters1["domain_model_part"].GetString())
solid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters1["problem_data"]["domain_size"].GetInt())

# import CFD_Thermal_solver           #we import the python file that includes the commands that we need
solver_module = __import__("CFD_thermal_solver")
solver_fluid = solver_module.CreateSolver(fluid_main_model_part,thermal_main_model_part, ProjectParameters)

solver_module1 = __import__("heat_equation_solver")
solver_solid = solver_module1.CreateSolver(solid_main_model_part, ProjectParameters1)


solver_fluid.AddVariables()
thermal_main_model_part.AddNodalSolutionStepVariable(FACE_HEAT_FLUX)
for node in thermal_main_model_part.Nodes:
    node.AddDof(TEMPERATURE,FACE_HEAT_FLUX)
fluid_main_model_part.AddNodalSolutionStepVariable(FACE_HEAT_FLUX)
for node in fluid_main_model_part.Nodes:
    node.AddDof(FACE_HEAT_FLUX)

solver_solid.AddVariables()

NonConformant_OneSideMap.AddVariables(fluid_main_model_part, thermal_main_model_part)

solver_fluid.ImportModelPart()
solver_solid.ImportModelPart()


 #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print rea results in a file
gid_mode = GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io_s = GidIO("art4_full_s",gid_mode,multifile,deformed_mesh_flag,write_conditions)
gid_io_f = GidIO("art4_full_f",gid_mode,multifile,deformed_mesh_flag,write_conditions)


 # we create a mesh for the postprocess

out_mesh_f = fluid_main_model_part.GetMesh()
out_mesh_s = solid_main_model_part.GetMesh()

out_nodes_f = fluid_main_model_part.Nodes
out_nodes_s = solid_main_model_part.Nodes

mesh_name_f = 0.0
mesh_name_s = 0.0
gid_io_f.InitializeMesh( mesh_name_f );
gid_io_s.InitializeMesh( mesh_name_s );
#
gid_io_f.WriteMesh( out_mesh_f );
gid_io_s.WriteMesh( out_mesh_s );
#
gid_io_f.FinalizeMesh()
gid_io_s.FinalizeMesh()


 # we add the DoFs

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transient problems, in this static problem =1 is enough)
fluid_main_model_part.SetBufferSize(3)
# print(fluid_main_model_part)
# thermal_main_model_part.SetBufferSize(3)
solid_main_model_part.SetBufferSize(3)
# print(solid_main_model_part)
# CalculateNodalAreaProcess(fluid_main_model_part,ProjectParameters["heat_settings"]["problem_data"]["domain_size"].GetInt()).Execute()


solver_fluid.AddDofs()
solver_solid.AddDofs()

solver_fluid.Initialize()
solver_solid.Initialize()



thermal_side = thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE")
solid_side = solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE")

fluid_nodes = thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Nodes
solid_nodes  = solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Nodes


fluid_nodes=list(fluid_nodes)
solid_nodes=list(solid_nodes)
#
id_list = ConnectivityMapper(fluid_nodes, solid_nodes)
print("number of matching nodes :", len(id_list))
interface_fluid_nodes = VariableUtils().SetFlag(INTERFACE, True, thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Nodes)
solid_interface_nodes = VariableUtils().SetFlag(INTERFACE, True, solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Nodes)

interface_fluid_condition = VariableUtils().SetFlag(INTERFACE, True, thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Conditions)
interface_solid_condition = VariableUtils().SetFlag(INTERFACE, True, solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Conditions)

interface_fluid_elements = VariableUtils().SetFlag(INTERFACE, True, thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Elements)
interface_solid_elements = VariableUtils().SetFlag(INTERFACE, True, solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Elements)

thermal_side = thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE")
solid_side = solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE")

# Set fluid and structure interfaces
for node in thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Nodes:
    node.Set(INTERFACE, True)
for node in solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Nodes:
    node.Set(INTERFACE, True)

# Mapper construction
search_radius_factor = 2.0
mapper_max_iterations = 200
mapper_tolerance = 1e-12
mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(thermal_side,
                                                                solid_side,
                                                                search_radius_factor,
                                                                mapper_max_iterations,mapper_tolerance)


# err
# Initial conditions
# Temperature fluid initial condition
initial_fluid_temperature = 298.15
#Temperature solid initial condition
initial_solid_temperature = 352.15

v_avarage = 2

    #fluid Temperature at the inlet
for node in thermal_main_model_part.GetSubModelPart("AutomaticInlet2D_INLET").Nodes:
    node.SetSolutionStepValue(TEMPERATURE,initial_fluid_temperature)
    node.Fix(TEMPERATURE)
    #fluid Temperature at t = 0
for node in thermal_main_model_part.Nodes:
    node.SetSolutionStepValue(TEMPERATURE,initial_fluid_temperature)

# outlet   V = 0, P = 0
for node in fluid_main_model_part.GetSubModelPart("Outlet2D_OUTLET").Nodes:
    Outlet_velocity = Vector(3)
    Outlet_velocity[1] = 0.0
    node.SetSolutionStepValue(VELOCITY, Outlet_velocity)
    node.Fix(VELOCITY_Y)
    node.SetSolutionStepValue(PRESSURE, 0.0)
    node.Fix(PRESSURE)

fluid_wall_velocity = Vector(3)
fluid_wall_velocity[0] = 0.0
fluid_wall_velocity[1] = 0.0
fluid_wall_velocity[2] = 0.0

# no-slip condition fluid
for node in fluid_main_model_part.GetSubModelPart("NoSlip2D_WALLS").Nodes:
    node.SetSolutionStepValue(VELOCITY, fluid_wall_velocity)
    node.Fix(VELOCITY_X)
    node.Fix(VELOCITY_Y)
    node.Fix(VELOCITY_Z)

# no-slip condition  fluid
for node in fluid_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Nodes:
    node.SetSolutionStepValue(VELOCITY, fluid_wall_velocity)
    node.Fix(VELOCITY_X)
    node.Fix(VELOCITY_Y)
    node.Fix(VELOCITY_Z)
# no-slip condition solid
for node in  solid_main_model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY, fluid_wall_velocity)
    node.Fix(VELOCITY_X)
    node.Fix(VELOCITY_Y)
    node.Fix(VELOCITY_Z)


for node in solid_main_model_part.Nodes:
    node.SetSolutionStepValue(HEAT_FLUX,0)
    node.Fix(HEAT_FLUX)

    # Temperature initial condition SOLID
for node in solid_main_model_part.Nodes:
    node.SetSolutionStepValue(TEMPERATURE,initial_solid_temperature)


time = 0.0
delta_time = 0.001
total_steps = 300000

only_fluid = 1/delta_time
all_model = only_fluid + 2

g0 = 9.81
gravity = Vector(3)
gravity[0] = 0.0
gravity[1] = g0
gravity[2] = 0.0

const = np.arcsin(1)/(only_fluid)*2/math.pi
for step in range(1,total_steps):

    time += delta_time
    print("Step: ", step)

    # if (step<=only_fluid):

    solid_main_model_part.CloneTimeStep(time)
    print("solid_main_model_part.CloneTimeStep(time) done!")
    solver_solid.SolverInitializeSolutionStep()
    print("solver_solid.SolverInitializeSolutionStep() OK")
    solver_solid.SolverPredict()
    print("solver_solid.SolverPredict() OK")
    if (step > 2):
        solver_solid.SolverSolveSolutionStep()
    solver_solid.SolverFinalizeSolutionStep()

    mapper.StructureToFluid_ScalarMap( TEMPERATURE, TEMPERATURE, True)

            #
            # k=0
            # fluid_therma_list = []
    for node in thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Nodes:
            # thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Nodes[id_list[k][0]].SetSolutionStepValue(TEMPERATURE, solid_temp_list[k]) # solid_temp_list[k]
        node.Fix(TEMPERATURE)
            #     k += 1


    fluid_main_model_part.CloneTimeStep(time) # solve the fluid
    print("fluid_main_model_part.CloneTimeStep(time) done!")
    for node in fluid_main_model_part.GetSubModelPart("AutomaticInlet2D_INLET").Nodes:
        Inlet_velocity = Vector(3)
        Inlet_velocity[0] = 4*v_avarage*node.Y*(0.2-node.Y)/(0.2*0.2)
        Inlet_velocity[1] = 0.0
        Inlet_velocity[2] = 0.0

        node.SetSolutionStepValue(VELOCITY, Inlet_velocity)
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        node.Fix(VELOCITY_Z)
    if (step > 2):
        print('before solve fluid')
        solver_fluid.Solve()
        print("solve fluid OK")
    # else:



        # i=0
        # solid_temp_list = []
        # for node in solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Nodes:
        #     solid_temp_list.append(solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Nodes[id_list[i][1]].GetSolutionStepValue(TEMPERATURE))
        #     i += 1
        # print('  TEMPERATURE  solution of the solid at the Interface')
        # print(solid_temp_list)





        # fluid_main_model_part.CloneTimeStep(time) # solve the fluid
        # print("fluid_main_model_part.CloneTimeStep(time) done!")
        #
        # for node in fluid_main_model_part.GetSubModelPart("AutomaticInlet2D_INLET").Nodes:
        #     Inlet_velocity = Vector(3)
        #     Inlet_velocity[0] = 4*v_avarage*node.Y*(0.2-node.Y)/(0.2*0.2)
        #     Inlet_velocity[1] = 0.0
        #     Inlet_velocity[2] = 0.0
        #     node.SetSolutionStepValue(VELOCITY, Inlet_velocity)
        #     node.Fix(VELOCITY_X)
        #     node.Fix(VELOCITY_Y)
        #     node.Fix(VELOCITY_Z)


    print('before solve fluid')
    for node in thermal_main_model_part.Nodes:
        Temperature = node.GetSolutionStepValue(TEMPERATURE)
        g = g0 * (1 - 1/initial_fluid_temperature * (Temperature - initial_fluid_temperature))
        # print(g)
        gravity[1] = g
        node.SetSolutionStepValue(BODY_FORCE, gravity)
        node.Fix(BODY_FORCE_X)
        node.Fix(BODY_FORCE_Y)
        node.Fix(BODY_FORCE_Z)
    solver_fluid.Solve() # solve the fluid
    print("solve fluid OK")
        #
    mapper.FluidToStructure_ScalarMap( FACE_HEAT_FLUX, FACE_HEAT_FLUX, True)
        #
        # fluid_flux_list = []
        # i=0
        # for node in thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Nodes:
        #     fluid_flux_list.append(thermal_main_model_part.GetSubModelPart("NoSlip2D_FLUID_INTERFACE").Nodes[id_list[i][0]].GetSolutionStepValue(FACE_HEAT_FLUX))
        #     i += 1
        # print(' fluid REACTION at the Interface')
        #
        # print(fluid_flux_list)
        # sleep(1)
        # #
        # k = 0
    for node in solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Nodes:
        #     solid_main_model_part.GetSubModelPart("NoSlip2D_SOLID_INTERFACE").Nodes[id_list[k][1]].SetSolutionStepValue(FACE_HEAT_FLUX,fluid_flux_list[k])
        node.Fix(FACE_HEAT_FLUX)
        #     k += 1


    print ("Solved!")
    gid_io_f.InitializeResults(time,out_mesh_f)
    gid_io_f.WriteNodalResults(BODY_FORCE,out_nodes_f,time,0)
    gid_io_f.WriteNodalResults(TEMPERATURE,out_nodes_f,time,0)
    gid_io_f.WriteNodalResults(VELOCITY,out_nodes_f,time,0)
    gid_io_f.WriteNodalResults(PRESSURE,out_nodes_f,time,0)
    gid_io_f.WriteNodalResults(FACE_HEAT_FLUX,out_nodes_f,time,0)
    gid_io_f.FinalizeResults()
    #
    gid_io_s.InitializeResults(time,out_mesh_s)
    gid_io_s.WriteNodalResults(TEMPERATURE,out_nodes_s,time,0)
    gid_io_s.WriteNodalResults(FACE_HEAT_FLUX,out_nodes_s,time,0)
    gid_io_s.WriteNodalResults(NODAL_AREA,out_nodes_s,time,0)
    gid_io_s.FinalizeResults()
