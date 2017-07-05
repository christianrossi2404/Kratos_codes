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



def AnalyticalSolution(node, a , time):
    return a * (node.X**2 + node.Y**2 + node.X*node.Y*time )

#defining a model part
model_part = ModelPart("ExampleModelPart");  #we create a model part

#ProjectParameters = Parameters("""{}""")
ProjectParameters = Parameters("""{
"solver_type": "heat_equation_solver",
"problem_data" :{
    "domain_size" : 2
    },
"solver_settings" : {
    "solver_type" : "heat_equation_solver",
    "model_import_settings": {
        "input_type": "mdpa",
        "input_filename": "SGS"
        },
    "time_order" : 2,
    "echo_level" : 0,
    "domain_model_part": "Parts_Domain",
     "dynamic_tau":1
    }
}""")
model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

settings = ProjectParameters["solver_settings"]

# import heat_equation_solver           #we import the python file that includes the commands that we need
# heat_equation_solver = heat_equation_solver.CreateSolver(model_part, settings)

solver_module = __import__(settings["solver_type"].GetString())
heat_equation_solver = solver_module.CreateSolver(model_part, settings)

heat_equation_solver.AddVariables()

# Compute nodal area to set source term
# model_part.AddNodalSolutionStepVariable(NODAL_AREA)

 # (note that our model part does not have nodes or elements yet)

 #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file
gid_mode = GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("art4",gid_mode,multifile,deformed_mesh_flag,write_conditions)


heat_equation_solver.ImportModelPart()

 # we create a mesh for the postprocess

out_mesh = model_part.GetMesh()
out_nodes = model_part.Nodes
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh(out_mesh);
gid_io.FinalizeMesh()

# gid_io.InitializeResults(0,out_mesh)
# gid_io.WriteNodalResults(HEAT_FLUX,out_nodes,0.0,0)
# gid_io.WriteNodalResults(NODAL_AREA,out_nodes,0.0,0)
# gid_io.WriteNodalResults(TEMPERATURE,out_nodes,0.0,0)
# gid_io.FinalizeResults()

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)
model_part.SetBufferSize(3)

 # we add the DoFs
heat_equation_solver.AddDofs()

heat_equation_solver.Initialize()


# CalculateNodalAreaProcess(model_part,ProjectParameters["problem_data"]["domain_size"].GetInt()).Execute()


conv_vel = Vector(3)
conv_vel[0] = 200
conv_vel[1] = 0.0
conv_vel[2] = 0.0
for node in model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY, conv_vel)

a = conv_vel[0]

# # BCs imposition (TODO: Use processes ASAP)

time_step = 0.1
L = 1.0/20.0

rho = 1
cp = 1
k = 1
alpha = k/(rho*cp)
Peclet = abs(a*L/(alpha))

time = 0.0
total_steps = 100

for step in range(1,total_steps):

    time += time_step
    model_part.CloneTimeStep(time)
    heat_equation_solver.SolverInitializeSolutionStep()

    heat_equation_solver.SolverPredict()

    for node in model_part.Nodes:
        heat_s = a*node.X*node.Y + a * (2*a*node.X + a*node.Y*time ) - (4*a)
        node.SetSolutionStepValue(HEAT_FLUX,heat_s)

    if (step >= 3):
        for node in model_part.GetSubModelPart("Inlet").Nodes:
             temp = AnalyticalSolution(node, a,time)
             node.Fix(TEMPERATURE)
             node.SetSolutionStepValue(TEMPERATURE, temp)
        for node in model_part.GetSubModelPart("Outlet").Nodes:
             temp = AnalyticalSolution(node, a,time)
             node.Fix(TEMPERATURE)
             node.SetSolutionStepValue(TEMPERATURE,temp)
        for node in model_part.GetSubModelPart("wall").Nodes:
             temp = AnalyticalSolution(node,a, time)
             node.Fix(TEMPERATURE)
             node.SetSolutionStepValue(TEMPERATURE, temp)

        print ("about to solve!")
        heat_equation_solver.SolverSolveSolutionStep()
        print ("Solved!")
    else:
    # Fill buffer with analytical solution
        for node in model_part.Nodes:
            temp = AnalyticalSolution(node, a,time)
            node.SetSolutionStepValue(TEMPERATURE, temp)

    heat_equation_solver.SolverFinalizeSolutionStep()

    #and we print the results
    # gid_io.InitializeResults(time,out_mesh)
    # gid_io.WriteNodalResults(HEAT_FLUX,out_nodes,time,0)
    # gid_io.WriteNodalResults(NODAL_AREA,out_nodes,time,0)
    # gid_io.WriteNodalResults(TEMPERATURE,out_nodes,time,0)
    # gid_io.FinalizeResults()


    gid_io.InitializeResults(time,out_mesh)
    gid_io.WriteNodalResults(HEAT_FLUX,out_nodes,time,0)
    gid_io.WriteNodalResults(NODAL_AREA,out_nodes,time,0)
    gid_io.WriteNodalResults(TEMPERATURE,out_nodes,time,0)
    gid_io.FinalizeResults()

    # Computational

h = np.arange(0.0,1.0+1.0/20.0,(1.0/20.0))
print(h.shape[0])

x = (h.shape[0])
time_step = np.arange(0.0,10,0.01)
# print(time_step)
t = (time_step.shape[0])
# print(t)
C = np.zeros((t,x))  # [10 x 20]

z = []
# Analytical
for i in range(1,t+1):  # Analyti solution of the points [X,0.45] for different time steps
    for j in range(1,x+1):
        z = a * (h[j-1]**2 + 0.45**2 + 0.45*h[j-1]*time_step[i-1])
        C[i-1][j-1] = z


# print(C[-1][:])

#Numerical
F = [40.5,
85.550003,
131.60001,
178.64999,
226.7,
275.75,
325.79999,
376.85001,
428.89999,
481.95001,
536,
591.04999,
647.09998,
704.15002,
762.20001,
821.25,
881.29999,
942.34998,
1004.4,
1067.45,
1131.5]

# shapeF = (F.shape[0])

plt.ylabel('Temperature [°C]')
plt.xlabel('Step size [m]')
    # plt.grid(True)

plt.plot(h,F,'blue',label="Numerical")

plt.plot(h,C[-1][:], 'ro',label="Analytical")

plt.legend(loc=2, borderaxespad=0.)
plt.show()
#
#
# CONDUCTIVITY 1 // 0.0257
# SPECIFIC_HEAT 1 // 1.005
# DENSITY 1 // 1.205
# print("CFD:" ,CFD)
print("Peclet:",Peclet  )
