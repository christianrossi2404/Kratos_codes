//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
//



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
//#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "heat_equation_application.h"
#include "includes/variables.h"


namespace Kratos
{
	//Example
// 	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

 	KratosHeatEquationApplication::KratosHeatEquationApplication():
		mHeatEquation2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
		mHeatEquation3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mHeatEquationNeumannCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mHeatEquationNeumannCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
	{}

 	void KratosHeatEquationApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
    KRATOS_REGISTER_ELEMENT("HeatEquation2D", mHeatEquation2D);
    KRATOS_REGISTER_ELEMENT("HeatEquation3D", mHeatEquation3D);
    KRATOS_REGISTER_CONDITION("HeatEquationNeumannCondition2D",mHeatEquationNeumannCondition2D);
    KRATOS_REGISTER_CONDITION("HeatEquationNeumannCondition3D",mHeatEquationNeumannCondition3D);
 		std::cout << "Initializing KratosHeatEquationApplication... " << std::endl;

// 		KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR )
// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);


 	}

}  // namespace Kratos.
