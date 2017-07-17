//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    Christian Rossi
#if !defined(KRATOS_CALCUALTE_COMPUTE_BOUSSINESQ_INCLUDED )
#define  KRATOS_CALCUALTE_COMPUTE_BOUSSINESQ_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "HeatEquationApplication.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	//template< unsigned int TDim>
	class ComputeBoussinesq: public Process
	{
	public:

		KRATOS_CLASS_POINTER_DEFINITION(ComputeBoussinesq);

		ComputeBoussinesq(ModelPart& model_part, Parameters& rParameters)
			: mr_model_part(model_part)              //mr_model_part is saved as private variable (declared at the end of the file)
		{

		}


		virtual~ComputeBoussinesq() override;
		{}


 		virtual void Execute() override;
		{
			// KRATOS_TRY

			// double area;                    //we create the needed variables
			// double sum_areas=0.0;
			// double sum_temperatures=0.0;
			// double one_third=1.0/3.0;

		// 	//getting data for the given geometry
		// 	for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); //looping the elements
		// 		ielem!=mr_model_part.ElementsEnd(); ielem++)
		// 	{
		// 		Geometry<Node<3> >& geom = ielem->GetGeometry();
		//
		// 		for (unsigned int k = 0; k < 3; k++)
		// 		{
		// 			double Temperature = geom[k].FastGetSolutionStepValue(TEMPERATURE)*one_third*area;
		// 		}
		// 	}
		//
		// 	KRATOS_CATCH("")
		// }
		};
		/// Turn back information as a string.
		virtual std::string Info() const override;

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const override;

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const override;


	protected:


	private:
		array_1d<double,3> mrGravity;

		ModelPart& mr_model_part;



}  // namespace Kratos.

#endif //   defined
