//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Christian Rossi
//

#include "custom_elements/heat_equation.h"

namespace Kratos {

  template<>
  void HeatEquation<3>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
  {
      KRATOS_TRY

      unsigned int Dim = 3;
      unsigned int NumNodes = 4;
      unsigned int DofSize  = NumNodes;

      if (rResult.size() != DofSize)
          rResult.resize(DofSize, false);

      for(unsigned int i=0; i<NumNodes; i++)
      {

          rResult[i]  =  this->GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
      }

      KRATOS_CATCH("")
  }


template<>
void HeatEquation<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int NumNodes = 3;
    unsigned int DofSize  = NumNodes;

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        rResult[i]  =  this->GetGeometry()[i].GetDof(TEMPERATURE).EquationId();

    }

    KRATOS_CATCH("")
}


template<>
void HeatEquation<3>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes;

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i ]  =  this->GetGeometry()[i].pGetDof(TEMPERATURE);

    }

    KRATOS_CATCH("");
}

template<>
void HeatEquation<2>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int NumNodes = 3;
    unsigned int DofSize  = NumNodes;

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {

        ElementalDofList[i]  =  this->GetGeometry()[i].pGetDof(TEMPERATURE);
    }

    KRATOS_CATCH("");
}




template<>
void HeatEquation<3>::ComputeGaussPointLHSContribution(bounded_matrix<double,4,4>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const double rho = data.rho;
    const double k = data.k;                              // thermal conductivity
    const double cp = data.cp;                                  //  specific heat at constant pressure

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    // const double& tau = data.tau;
    const double h = data.h;                                // Characteristic element size
    const double& delta_t = data.delta_t;
    const array_1d<double,nnodes>& Q = data.Q;                              //Source term
    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const array_1d<double,nnodes>& temp = data.temp;
    const array_1d<double,nnodes>& tempn = data.tempn;
    const array_1d<double,nnodes>& tempnn = data.tempnn;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const double v_norm = norm_2(v_gauss);
    // Stabilization parameters
    const double c1 = 2.0;
    const double c2 = 4.0;
    const double tau = 1.0/(dyn_tau_coeff*rho*cp/delta_t + c1*rho*cp*v_norm/h + c2*k/(h*h));
    // const double tau = dyn_tau_coeff*1.0/(c1*rho*cp*v_norm/h + c2*k/(h*h));
    //substitute_lhs_3D

}


template<>
void HeatEquation<2>::ComputeGaussPointLHSContribution(bounded_matrix<double,3,3>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
	  const double rho = data.rho;
	  const double k = data.k;
    const double cp = data.cp;
   const double h = data.h;                                // Characteristic element size
   const double& delta_t = data.delta_t;
   const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    // const double& tau = data.tau;
    const array_1d<double,nnodes>& Q = data.Q;   //Source term
    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const array_1d<double,nnodes>& temp = data.temp;
    const array_1d<double,nnodes>& tempn = data.tempn;
    const array_1d<double,nnodes>& tempnn = data.tempnn;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;



    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    // KRATOS_WATCH(N)
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
    // KRATOS_WATCH(DN)

    // Stabilization parameters
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const double v_norm = norm_2(v_gauss);
    // Stabilization parameters
    const double c1 = 2.0;
    const double c2 = 4.0;
    const double tau = 1.0/(dyn_tau_coeff*rho*cp/delta_t + c1*rho*cp*v_norm/h + c2*k/(h*h));


}


template<>
void HeatEquation<3>::ComputeGaussPointRHSContribution(array_1d<double,4>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const double rho = data.rho;
    const double k = data.k;                              // thermal conductivity
    const double cp = data.cp;                            //  specific heat at constant pressure
    const double h = data.h;                                // Characteristic element size
    const double& delta_t = data.delta_t;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    //const double& tau = data.tau;
    const array_1d<double,nnodes>& Q = data.Q;   //Source term
    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const array_1d<double,nnodes>& temp = data.temp;
    const array_1d<double,nnodes>& tempn = data.tempn;
    const array_1d<double,nnodes>& tempnn = data.tempnn;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const double q_gauss = inner_prod(data.Q, data.N);

    //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
    //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
    //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);

    // Stabilization parameters
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const double v_norm = norm_2(v_gauss);
    // Stabilization parameters
    const double c1 = 2.0;
    const double c2 = 4.0;
    const double tau = 1.0/(dyn_tau_coeff*rho*cp/delta_t + c1*rho*cp*v_norm/h + c2*k/(h*h));
    // const double tau = dyn_tau_coeff*1.0/(c1*rho*cp*v_norm/h + c2*k/(h*h));

    //substitute_rhs_3D
}


template<>
void HeatEquation<2>::ComputeGaussPointRHSContribution(array_1d<double,3>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    const double rho = data.rho;
    const double k = data.k;                              // thermal conductivity
    const double cp = data.cp;                            //  specific heat at constant pressure
    const double h = data.h;                                // Characteristic element size
    const double& delta_t = data.delta_t;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    //const double& tau = data.tau;
    const array_1d<double,nnodes>& Q = data.Q;  //Source term
    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const array_1d<double,nnodes>& temp = data.temp;
    const array_1d<double,nnodes>& tempn = data.tempn;
    const array_1d<double,nnodes>& tempnn = data.tempnn;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS

    const double q_gauss = inner_prod(data.Q, data.N);
    //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
    //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
    //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);

    // Stabilization parameters
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const double v_norm = norm_2(v_gauss);
    // Stabilization parameters
    const double c1 = 2.0;
    const double c2 = 4.0;
    const double tau = 1.0/(dyn_tau_coeff*rho*cp/delta_t + c1*rho*cp*v_norm/h + c2*k/(h*h));
    // const double tau = dyn_tau_coeff*1.0/(c1*rho*cp*v_norm/h + c2*k/(h*h));
    //substitute_rhs_2D
    // KRATOS_WATCH(tau)
}

}
