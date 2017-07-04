from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.HeatEquationApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(fluid_model_part, solid_model_part):
    fluid_model_part.AddNodalSolutionStepVariable(NODAL_MAUX)  # Stores Nodal Area
    fluid_model_part.AddNodalSolutionStepVariable(MAPPER_SCALAR_PROJECTION_RHS)
    fluid_model_part.AddNodalSolutionStepVariable(MAPPER_VECTOR_PROJECTION_RHS)
    fluid_model_part.AddNodalSolutionStepVariable(VAUX_EQ_TRACTION)
    fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    fluid_model_part.AddNodalSolutionStepVariable(NORMAL)
    fluid_model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED)


    solid_model_part.AddNodalSolutionStepVariable(NODAL_MAUX)  # Stores Nodal Area
    solid_model_part.AddNodalSolutionStepVariable(PRESSURE)
    solid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    solid_model_part.AddNodalSolutionStepVariable(MAPPER_SCALAR_PROJECTION_RHS)
    solid_model_part.AddNodalSolutionStepVariable(MAPPER_VECTOR_PROJECTION_RHS)
    solid_model_part.AddNodalSolutionStepVariable(VAUX_EQ_TRACTION)
    solid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    solid_model_part.AddNodalSolutionStepVariable(NORMAL)
    solid_model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED)


    print("Mapper variables added correctly.")


class NonConformant_OneSideMap:

    def __init__(self, fluid_model_part, solid_model_part,
                 search_radius_factor=2.0, it_max=25, tol=1e-5):

        # Error check
        # if fluid_model_part.NumberOfConditions(0) < 1:
        #     raise ValueError("No conditions found in the fluid model part, please check that the interface is meshed using SurfaceCondition2D2N/SurfaceCondition3D3N")
        # if solid_model_part.NumberOfConditions(0) < 1:
        #     raise ValueError("No conditions found in the solid model part, please check that the interface is meshed using SurfaceCondition2D2N/SurfaceCondition3D3N")

        search_radius_factor = search_radius_factor
        self.it_max = it_max
        self.tol = tol

        self.Preprocess = InterfacePreprocess()
        self.fl_interface = ModelPart("fluid_interface")
        self.str_interface = ModelPart("solid_interface")

        domain_size_fl = fluid_model_part.ProcessInfo[DOMAIN_SIZE]
        domain_size_solid = solid_model_part.ProcessInfo[DOMAIN_SIZE]

        if (domain_size_fl != domain_size_fl):
          raise ValueError("Domain sizes from two model parts are not compatible")

        print("Identifying fluid interface...")
        if (domain_size_fl == 3):
          self.Preprocess.GenerateTriangleInterfacePart(fluid_model_part, self.fl_interface)
        else:
          self.Preprocess.GenerateLineInterfacePart(fluid_model_part, self.fl_interface)

        print("Identifying solid interface...")
        if (domain_size_fl == 3):
          self.Preprocess.GenerateTriangleInterfacePart(solid_model_part, self.str_interface)
        else:
          self.Preprocess.GenerateLineInterfacePart(solid_model_part, self.str_interface)
        print("Fluid and solid interfaces identified.")

        # Interface mappers construction
        self.FluidToStructureSolid = AdvancedNMPointsMapper\
            (self.fl_interface, self.str_interface)
        self.SolidToFluidMapper = AdvancedNMPointsMapper\
            (self.str_interface, self.fl_interface)
        print("Interface Mappers created.")

        # Neighbour search
        (self.FluidToStructureSolid).FindNeighbours(search_radius_factor)
        (self.SolidToFluidMapper).FindNeighbours(search_radius_factor)
        print("Neighbours search finished.")

    def RecomputeTransferPairs(self, search_radius_factor):
        (self.FluidToStructureSolid).FindNeighbours(search_radius_factor)
        (self.SolidToFluidMapper).FindNeighbours(search_radius_factor)

    # Standard mappers
    def SolidToFluid_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        (self.SolidToFluidMapper).VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos, distributed)

    def SolidToFluid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.SolidToFluidMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def FluidToSolid_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        (self.FluidToStructureSolid).VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos, distributed)

    def FluidToSolid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.FluidToStructureSolid).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    # Normal vectors
    def StructureToFluid_ScalarToNormalVectorMap(self, ScalarVar_Origin, VectorVar_Destination, sign_pos):
        (self.SolidToFluidMapper).ScalarToNormalVectorMap(ScalarVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def StructureToFluid_NormalVectorToScalarMap(self, VectorVar_Origin, ScalarVar_Destination, sign_pos):
        (self.SolidToFluidMapper).NormalVectorToScalarMap(VectorVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def FluidToStructure_ScalarToNormalVectorMap(self, ScalarVar_Origin, VectorVar_Destination, sign_pos):
        (self.FluidToStructureSolid).ScalarToNormalVectorMap(ScalarVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def FluidToStructure_NormalVectorToScalarMap(self, VectorVar_Origin, ScalarVar_Destination, sign_pos):
        (self.FluidToStructureSolid).NormalVectorToScalarMap(VectorVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)
