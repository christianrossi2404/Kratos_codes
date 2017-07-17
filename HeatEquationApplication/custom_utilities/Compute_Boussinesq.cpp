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


#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "custom_utilities/ComputeBoussinesq.h"


namespace Kratos
{




  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;


		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

		class_<ComputeBoussinesq > ("ComputeBoussinesq", init<ModelPart& >())  //the input parameters is a model part
                   .def("Execute", &ComputeBoussinesq::Calculate)  //when we call "Execute" in python, Calculate is called in C++. Notice we don't write the input parameters here
                   ;


  }





} // Namespace Kratos
