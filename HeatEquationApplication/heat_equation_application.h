//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_HEATEQUATION_APPLICATION_H_INCLUDED )
#define  KRATOS_HEATEQUATION_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "custom_conditions/heat_equation_Neumann_condition.h"
#include "custom_elements/heat_equation.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

	///@name Kratos Globals
	///@{

	// Variables definition
//	KRATOS_DEFINE_VARIABLE(double, AUX_MESH_VAR )
//	KRATOS_DEFINE_VARIABLE(double, IS_INTERFACE)
//	KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)


	///@}
	///@name Type Definitions
	///@{

	///@}
	///@name  Enum's
	///@{

	///@}
	///@name  Functions
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
	*/
	class KratosHeatEquationApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{


		/// Pointer definition of KratosHeatEquationApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosHeatEquationApplication);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		KratosHeatEquationApplication();

		/// Destructor.
		virtual ~KratosHeatEquationApplication(){}


		///@}
		///@name Operators
		///@{


		///@}
		///@name Operations
		///@{

		virtual void Register();



		///@}
		///@name Access
		///@{


		///@}
		///@name Inquiry
		///@{


		///@}
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosHeatEquationApplication";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
      	KRATOS_WATCH("in my application");
      	KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
      }


		///@}
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables
		///@{


		///@}
		///@name Protected member Variables
		///@{


		///@}
		///@name Protected Operators
		///@{


		///@}
		///@name Protected Operations
		///@{


		///@}
		///@name Protected  Access
		///@{


		///@}
		///@name Protected Inquiry
		///@{


		///@}
		///@name Protected LifeCycle
		///@{


		///@}

	private:
		///@name Static Member Variables
		///@{

		const HeatEquation<2>  mHeatEquation2D;
		const HeatEquation<3>  mHeatEquation3D;
		const HeatEquationNeumannCondition<2> mHeatEquationNeumannCondition2D;
		const HeatEquationNeumannCondition<3> mHeatEquationNeumannCondition3D;
		//       static const ApplicationCondition  msApplicationCondition;

		///@}
		///@name Member Variables
		///@{
// 		const Elem2D   mElem2D;
// 		const Elem3D   mElem3D;


		///@}
		///@name Private Operators
		///@{


		///@}
		///@name Private Operations
		///@{


		///@}
		///@name Private  Access
		///@{


		///@}
		///@name Private Inquiry
		///@{


		///@}
		///@name Un accessible methods
		///@{

		/// Assignment operator.
		KratosHeatEquationApplication& operator=(KratosHeatEquationApplication const& rOther);

		/// Copy constructor.
		KratosHeatEquationApplication(KratosHeatEquationApplication const& rOther);


		///@}

	}; // Class KratosHeatEquationApplication

	///@}


	///@name Type Definitions
	///@{


	///@}
	///@name Input and output
	///@{

	///@}


}  // namespace Kratos.

#endif // KRATOS_HEATEQUATION_APPLICATION_H_INCLUDED  defined
