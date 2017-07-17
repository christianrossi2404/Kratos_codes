//
//  Main authors:    Miguel Ángel Celigueta
//
//

#ifndef KRATOS_COMPUTE_BODY_FORCE_PROCESS_H
#define KRATOS_COMPUTE_BODY_FORCE_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"

// Application includes
namespace Kratos
{

class ComputeBodyForceProcess : public Process
{
public:

    /// Pointer definition of ComputeBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeBodyForceProcess);

    typedef Node<3>                     NodeType;
    typedef Geometry<NodeType>      GeometryType;

    /// Constructor.
    ComputeBodyForceProcess(ModelPart& rModelPart,
                    const double initial_fluid_temperature
                    const double gravity_x
                    const double gravity_y
                    const double gravity_z
                    const double fluid_thermal_expansion
    {
        mGravity[0] = 0.0;
        mGravity[1] = gravity_y;
        mGravity[2] = 0.0;
        mInitial_fluid_temperature = initial_fluid_temperature;
        mFluid_thermal_expansion = 1/initial_fluid_temperature;
    }

    ComputeBodyForceProcess(ModelPart& rModelPart,
                    Parameters rParameters): Process(), mrModelPart(rModelPart)
    {
        Parameters default_parameters( R"(
                {
                    "model_part_name"                                   :"model_part",
                    "initial_fluid_temperature"                         :298.15,
                    "gravity_x"                                         :0.0,
                    "gravity_y"                                         :9.81,
                    "gravity_z"                                         :0.0,
                    "fluid_thermal_expansion"                           :""
                } )"   );

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mGravity[0] = 0.0;
        mGravity[1] = rParameters["gravity_y"].GetDouble();
        mGravity[2] = 0.0;
        mInitial_fluid_temperature = rParameters["initial_fluid_temperature"].GetDouble();
        mFluid_thermal_expansion = 1/initial_fluid_temperature;

    }

    /// Destructor.
    virtual ~ComputeBodyForceProcess(){}


    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_TRY;



        KRATOS_CATCH("");
    }


    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
	    const double& rCurrentTime = rCurrentProcessInfo[TIME];

        if (!mrModelPart.NodesBegin()->SolutionStepsDataHas(TEMPERATURE
            KRATOS_ERROR << "TEMPERATURE variable is not in submodelpart.";

        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_1, vertical_vector, rotated_angle2, current_local_axis_1);
        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_2, vertical_vector, rotated_angle2, current_local_axis_2);


        int NumNodes = rModelPart.NumberOfNodes();
        #pragma omp parallel for firstprivate(NumNodes,initial_fluid_temperature)
        for (int i = 0; i < NumNodes; ++i)
        {
            ModelPart::NodeIterator iNode = rModelPart.NodesBegin() + i;
            double Temperature = iNode->FastGetSolutionStepValue(TEMPERATURE);

            iNode->FastGetSolutionStepValue(BODY_FORCE) = (1. - fluid_thermal_expansion *(Temperature-initial_fluid_temperature))*mrGravity;
        }
        //UPDATE POSITION AND VELOCITY OF ALL NODES
        for (ModelPart::NodesContainerType::iterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); node_i++) {
            //Get local coordinates at the beginning (local axes are assumed oriented as global axes at the beginning)
            array_1d<double,3> local_coordinates;
            local_coordinates[0] = node_i->X0() - mInitialCoordinatesOfRotorCenter[0];
            local_coordinates[1] = node_i->Y0() - mInitialCoordinatesOfRotorCenter[1];

            //Use local coordinates with the updated local axes
            array_1d<double,3> from_rotor_center_to_node;
            noalias(from_rotor_center_to_node) = local_coordinates[0] * current_local_axis_1 + local_coordinates[1] * current_local_axis_2;

            array_1d<double,3>& current_node_position = node_i->Coordinates();
            current_node_position[0] = current_rotor_position[0] + local_coordinates[0] * current_local_axis_1[0] + local_coordinates[1] * current_local_axis_2[0];
            current_node_position[1] = current_rotor_position[1] + local_coordinates[0] * current_local_axis_1[1] + local_coordinates[1] * current_local_axis_2[1];

            noalias(node_i->FastGetSolutionStepValue(DISPLACEMENT)) = node_i->Coordinates() - node_i->GetInitialPosition();

            array_1d<double,3>& current_node_velocity = node_i->FastGetSolutionStepValue(VELOCITY);
            noalias(current_node_velocity) = rotor_velocity + MathUtils<double>::CrossProduct(from_rotor_center_to_node, mW2);

        }//end of loop over nodes
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ComputeBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ComputeBodyForceProcess";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart&                               mrModelPart;
    array_1d<double,3>                       mGravity;

    double                                   mInitial_fluid_temperature;
    double                                   mFluid_thermal_expansion;

    }



private:

    /// Assignment operator.
    ComputeBodyForceProcess& operator=(ComputeBodyForceProcess const& rOther){return *this;}

}; // Class ComputeBodyForceProcess

};  // namespace Kratos.

#endif // KRATOS_COMPUTE_BODY_FORCE_PROCESS_H
