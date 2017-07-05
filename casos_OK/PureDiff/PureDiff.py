##########################################################
##################################################################
#setting the domain size for the problem to be solved

#including kratos path
import sys
import numpy as np
import matplotlib.pyplot as plt
from KratosMultiphysics import *    #we import the KRATOS
from KratosMultiphysics.HeatEquationApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem
import KratosMultiphysics.ExternalSolversApplication
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
        "input_filename": "PureDiff"
        },
    "time_order" : 2,
    "echo_level" : 0,
    "domain_model_part": "Parts_Domain"
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
# CalculateNodalAreaProcess(model_part,ProjectParameters["problem_data"]["domain_size"].GetInt()).Execute()

 # (note that our model part does not have nodes or elements yet)

 #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file
gid_mode = GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("art4",gid_mode,multifile,deformed_mesh_flag,write_conditions)

model_part_io = ModelPartIO("PureDiff")             # we set the name of the .mdpa file
model_part_io.ReadModelPart(model_part)         # we load the info from the .mdpa

 # we create a mesh for the postprocess

out_mesh = model_part.GetMesh()
out_nodes = model_part.Nodes

gid_io.InitializeMesh( 0.0 );
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
# for node in model_part.Nodes:
#     node.SetSolutionStepValue(HEAT_FLUX, -4.0)

time = 0.0

for step in range(1,10):

    time += 0.1
    model_part.CloneTimeStep(time)

    for node in model_part.Nodes:
        heat_s = node.X*node.Y  - 4.0
        node.SetSolutionStepValue(HEAT_FLUX,heat_s)

    # BCs imposition (TODO: Use processes ASAP)
    for node in model_part.GetSubModelPart("AutomaticInlet2D_Boundary").Nodes:
        temp = node.X**2 + node.Y**2 + node.X*node.Y*time
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE, temp)

    heat_equation_solver.SolverInitializeSolutionStep()

    heat_equation_solver.SolverPredict()

    if (step >= 3):
        print ("about to solve!")
        heat_equation_solver.SolverSolveSolutionStep()
        print ("Solved!")
    else:
        # Fill buffer with analytical solution
        for node in model_part.Nodes:
            temp = node.X**2 + node.Y**2 + node.X*node.Y*time
            node.SetSolutionStepValue(TEMPERATURE, temp)

    heat_equation_solver.SolverFinalizeSolutionStep()

    #and we print the results
    gid_io.InitializeResults(time,out_mesh)
    gid_io.WriteNodalResults(HEAT_FLUX,out_nodes,time,0)
    gid_io.WriteNodalResults(NODAL_AREA,out_nodes,time,0)
    gid_io.WriteNodalResults(TEMPERATURE,out_nodes,time,0)
    gid_io.FinalizeResults()

h = np.arange(0.0,1.0+1.0/3.0,1.0/3.0) # space discretization
time_step = np.arange(0.0,1.0,0.1) # time discretization
x = (h.shape[0])
t = (time_step.shape[0])
C = np.zeros((t,x))  # [10 x 20]
z = []
# Analytical
yy = 0.333333
for i in range(1,t+1):  # Analyti solution of the points [X,0.45] for different time steps
    for j in range(1,x+1):
        z = h[j-1]**2 + yy**2 + yy*h[j-1]*time_step[i-1]
        C[i-1][j-1] = z

# plt.title('Pure Diffusion at 0.9 s')
plt.ylabel('Temperature [°C]')
plt.xlabel('Step size')
F = [0.111111,0.3222222,0.7555556,1.411111]

# plt.grid(True)
plt.plot(h,C[9][:], 'ro',label="Analytical")
plt.plot(h,F,label="Numerical")
plt.legend(loc=2, borderaxespad=0.)
plt.show()
