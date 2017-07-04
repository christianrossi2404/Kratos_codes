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
    const double tau = dyn_tau_coeff*1.0/(rho*cp/delta_t + c1*rho*cp*v_norm/h + c2*k/(h*h));
    // const double tau = dyn_tau_coeff*1.0/(c1*rho*cp*v_norm/h + c2*k/(h*h));
    const double clhs0 =             bdf0*cp*rho;
const double clhs1 =             N[0]*cp*rho;
const double clhs2 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double clhs3 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double clhs4 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double clhs5 =             DN(0,0)*clhs2 + DN(0,1)*clhs3 + DN(0,2)*clhs4;
const double clhs6 =             N[0]*bdf0 + clhs5;
const double clhs7 =             pow(cp, 2);
const double clhs8 =             pow(rho, 2);
const double clhs9 =             clhs5*clhs7*clhs8*tau;
const double clhs10 =             DN(0,0)*k;
const double clhs11 =             DN(0,1)*k;
const double clhs12 =             DN(0,2)*k;
const double clhs13 =             N[0]*bdf0*cp*rho;
const double clhs14 =             DN(1,0)*clhs10 + DN(1,1)*clhs11 + DN(1,2)*clhs12 + N[1]*clhs13;
const double clhs15 =             DN(1,0)*clhs2 + DN(1,1)*clhs3 + DN(1,2)*clhs4;
const double clhs16 =             N[1]*bdf0 + clhs15;
const double clhs17 =             DN(2,0)*clhs10 + DN(2,1)*clhs11 + DN(2,2)*clhs12 + N[2]*clhs13;
const double clhs18 =             DN(2,0)*clhs2 + DN(2,1)*clhs3 + DN(2,2)*clhs4;
const double clhs19 =             N[2]*bdf0;
const double clhs20 =             clhs18 + clhs19;
const double clhs21 =             DN(3,0)*clhs10 + DN(3,1)*clhs11 + DN(3,2)*clhs12 + N[3]*clhs13;
const double clhs22 =             DN(3,0)*clhs2 + DN(3,1)*clhs3 + DN(3,2)*clhs4;
const double clhs23 =             N[3]*bdf0 + clhs22;
const double clhs24 =             N[1]*cp*rho;
const double clhs25 =             clhs15*clhs7*clhs8*tau;
const double clhs26 =             DN(1,0)*k;
const double clhs27 =             DN(1,1)*k;
const double clhs28 =             DN(1,2)*k;
const double clhs29 =             N[1]*bdf0*cp*rho;
const double clhs30 =             DN(2,0)*clhs26 + DN(2,1)*clhs27 + DN(2,2)*clhs28 + N[2]*clhs29;
const double clhs31 =             DN(3,0)*clhs26 + DN(3,1)*clhs27 + DN(3,2)*clhs28 + N[3]*clhs29;
const double clhs32 =             N[2]*cp*rho;
const double clhs33 =             clhs18*clhs7*clhs8*tau;
const double clhs34 =             N[3]*cp*rho;
const double clhs35 =             DN(2,0)*DN(3,0)*k + DN(2,1)*DN(3,1)*k + DN(2,2)*DN(3,2)*k + clhs19*clhs34;
const double clhs36 =             clhs22*clhs7*clhs8*tau;
            lhs(0,0)=pow(DN(0,0), 2)*k + pow(DN(0,1), 2)*k + pow(DN(0,2), 2)*k + pow(N[0], 2)*clhs0 + clhs1*clhs5 + clhs6*clhs9;
            lhs(0,1)=clhs1*clhs15 + clhs14 + clhs16*clhs9;
            lhs(0,2)=clhs1*clhs18 + clhs17 + clhs20*clhs9;
            lhs(0,3)=clhs1*clhs22 + clhs21 + clhs23*clhs9;
            lhs(1,0)=clhs14 + clhs24*clhs5 + clhs25*clhs6;
            lhs(1,1)=pow(DN(1,0), 2)*k + pow(DN(1,1), 2)*k + pow(DN(1,2), 2)*k + pow(N[1], 2)*clhs0 + clhs15*clhs24 + clhs16*clhs25;
            lhs(1,2)=clhs18*clhs24 + clhs20*clhs25 + clhs30;
            lhs(1,3)=clhs22*clhs24 + clhs23*clhs25 + clhs31;
            lhs(2,0)=clhs17 + clhs32*clhs5 + clhs33*clhs6;
            lhs(2,1)=clhs15*clhs32 + clhs16*clhs33 + clhs30;
            lhs(2,2)=pow(DN(2,0), 2)*k + pow(DN(2,1), 2)*k + pow(DN(2,2), 2)*k + pow(N[2], 2)*clhs0 + clhs18*clhs32 + clhs20*clhs33;
            lhs(2,3)=clhs22*clhs32 + clhs23*clhs33 + clhs35;
            lhs(3,0)=clhs21 + clhs34*clhs5 + clhs36*clhs6;
            lhs(3,1)=clhs15*clhs34 + clhs16*clhs36 + clhs31;
            lhs(3,2)=clhs18*clhs34 + clhs20*clhs36 + clhs35;
            lhs(3,3)=pow(DN(3,0), 2)*k + pow(DN(3,1), 2)*k + pow(DN(3,2), 2)*k + pow(N[3], 2)*clhs0 + clhs22*clhs34 + clhs23*clhs36;


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
    const double tau = dyn_tau_coeff*1.0/(rho*cp/delta_t + c1*rho*cp*v_norm/h + c2*k/(h*h));
    // const double tau = dyn_tau_coeff*1.0/(c1*rho*cp*v_norm/h + c2*k/(h*h));
    // KRATOS_WATCH(tau)
    const double clhs0 =             bdf0*cp*rho;
const double clhs1 =             N[0]*cp*rho;
const double clhs2 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double clhs3 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double clhs4 =             DN(0,0)*clhs2 + DN(0,1)*clhs3;
const double clhs5 =             N[0]*bdf0 + clhs4;
const double clhs6 =             pow(cp, 2);
const double clhs7 =             pow(rho, 2);
const double clhs8 =             clhs4*clhs6*clhs7*tau;
const double clhs9 =             DN(0,0)*k;
const double clhs10 =             DN(0,1)*k;
const double clhs11 =             N[0]*bdf0*cp*rho;
const double clhs12 =             DN(1,0)*clhs9 + DN(1,1)*clhs10 + N[1]*clhs11;
const double clhs13 =             DN(1,0)*clhs2 + DN(1,1)*clhs3;
const double clhs14 =             N[1]*bdf0;
const double clhs15 =             clhs13 + clhs14;
const double clhs16 =             DN(2,0)*clhs9 + DN(2,1)*clhs10 + N[2]*clhs11;
const double clhs17 =             DN(2,0)*clhs2 + DN(2,1)*clhs3;
const double clhs18 =             N[2]*bdf0 + clhs17;
const double clhs19 =             N[1]*cp*rho;
const double clhs20 =             clhs13*clhs6*clhs7*tau;
const double clhs21 =             N[2]*cp*rho;
const double clhs22 =             DN(1,0)*DN(2,0)*k + DN(1,1)*DN(2,1)*k + clhs14*clhs21;
const double clhs23 =             clhs17*clhs6*clhs7*tau;
            lhs(0,0)=pow(DN(0,0), 2)*k + pow(DN(0,1), 2)*k + pow(N[0], 2)*clhs0 + clhs1*clhs4 + clhs5*clhs8;
            lhs(0,1)=clhs1*clhs13 + clhs12 + clhs15*clhs8;
            lhs(0,2)=clhs1*clhs17 + clhs16 + clhs18*clhs8;
            lhs(1,0)=clhs12 + clhs19*clhs4 + clhs20*clhs5;
            lhs(1,1)=pow(DN(1,0), 2)*k + pow(DN(1,1), 2)*k + pow(N[1], 2)*clhs0 + clhs13*clhs19 + clhs15*clhs20;
            lhs(1,2)=clhs17*clhs19 + clhs18*clhs20 + clhs22;
            lhs(2,0)=clhs16 + clhs21*clhs4 + clhs23*clhs5;
            lhs(2,1)=clhs13*clhs21 + clhs15*clhs23 + clhs22;
            lhs(2,2)=pow(DN(2,0), 2)*k + pow(DN(2,1), 2)*k + pow(N[2], 2)*clhs0 + clhs17*clhs21 + clhs18*clhs23;


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
    const double tau = dyn_tau_coeff*1.0/(rho*cp/delta_t + c1*rho*cp*v_norm/h + c2*k/(h*h));
    // const double tau = dyn_tau_coeff*1.0/(c1*rho*cp*v_norm/h + c2*k/(h*h));

    const double crhs0 =             N[0]*Q[0] + N[1]*Q[1] + N[2]*Q[2] + N[3]*Q[3];
const double crhs1 =             DN(0,0)*temp[0] + DN(1,0)*temp[1] + DN(2,0)*temp[2] + DN(3,0)*temp[3];
const double crhs2 =             crhs1*k;
const double crhs3 =             DN(0,1)*temp[0] + DN(1,1)*temp[1] + DN(2,1)*temp[2] + DN(3,1)*temp[3];
const double crhs4 =             crhs3*k;
const double crhs5 =             DN(0,2)*temp[0] + DN(1,2)*temp[1] + DN(2,2)*temp[2] + DN(3,2)*temp[3];
const double crhs6 =             crhs5*k;
const double crhs7 =             N[0]*(bdf0*temp[0] + bdf1*tempn[0] + bdf2*tempnn[0]) + N[1]*(bdf0*temp[1] + bdf1*tempn[1] + bdf2*tempnn[1]) + N[2]*(bdf0*temp[2] + bdf1*tempn[2] + bdf2*tempnn[2]) + N[3]*(bdf0*temp[3] + bdf1*tempn[3] + bdf2*tempnn[3]);
const double crhs8 =             cp*crhs7*rho;
const double crhs9 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs10 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs11 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs12 =             crhs1*crhs9 + crhs10*crhs3 + crhs11*crhs5;
const double crhs13 =             cp*crhs12*rho;
const double crhs14 =             cp*rho*tau*(-cp*rho*(crhs12 + crhs7) + crhs0);
            rhs[0]=-DN(0,0)*crhs2 - DN(0,1)*crhs4 - DN(0,2)*crhs6 + N[0]*crhs0 - N[0]*crhs13 - N[0]*crhs8 + crhs14*(DN(0,0)*crhs9 + DN(0,1)*crhs10 + DN(0,2)*crhs11);
            rhs[1]=-DN(1,0)*crhs2 - DN(1,1)*crhs4 - DN(1,2)*crhs6 + N[1]*crhs0 - N[1]*crhs13 - N[1]*crhs8 + crhs14*(DN(1,0)*crhs9 + DN(1,1)*crhs10 + DN(1,2)*crhs11);
            rhs[2]=-DN(2,0)*crhs2 - DN(2,1)*crhs4 - DN(2,2)*crhs6 + N[2]*crhs0 - N[2]*crhs13 - N[2]*crhs8 + crhs14*(DN(2,0)*crhs9 + DN(2,1)*crhs10 + DN(2,2)*crhs11);
            rhs[3]=-DN(3,0)*crhs2 - DN(3,1)*crhs4 - DN(3,2)*crhs6 + N[3]*crhs0 - N[3]*crhs13 - N[3]*crhs8 + crhs14*(DN(3,0)*crhs9 + DN(3,1)*crhs10 + DN(3,2)*crhs11);

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
    const double tau = dyn_tau_coeff*1.0/(rho*cp/delta_t + c1*rho*cp*v_norm/h + c2*k/(h*h));
    // const double tau = dyn_tau_coeff*1.0/(c1*rho*cp*v_norm/h + c2*k/(h*h));
    const double crhs0 =             N[0]*Q[0] + N[1]*Q[1] + N[2]*Q[2];
const double crhs1 =             DN(0,0)*temp[0] + DN(1,0)*temp[1] + DN(2,0)*temp[2];
const double crhs2 =             crhs1*k;
const double crhs3 =             DN(0,1)*temp[0] + DN(1,1)*temp[1] + DN(2,1)*temp[2];
const double crhs4 =             crhs3*k;
const double crhs5 =             N[0]*(bdf0*temp[0] + bdf1*tempn[0] + bdf2*tempnn[0]) + N[1]*(bdf0*temp[1] + bdf1*tempn[1] + bdf2*tempnn[1]) + N[2]*(bdf0*temp[2] + bdf1*tempn[2] + bdf2*tempnn[2]);
const double crhs6 =             cp*crhs5*rho;
const double crhs7 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crhs8 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs9 =             crhs1*crhs7 + crhs3*crhs8;
const double crhs10 =             cp*crhs9*rho;
const double crhs11 =             cp*rho*tau*(-cp*rho*(crhs5 + crhs9) + crhs0);
            rhs[0]=-DN(0,0)*crhs2 - DN(0,1)*crhs4 + N[0]*crhs0 - N[0]*crhs10 - N[0]*crhs6 + crhs11*(DN(0,0)*crhs7 + DN(0,1)*crhs8);
            rhs[1]=-DN(1,0)*crhs2 - DN(1,1)*crhs4 + N[1]*crhs0 - N[1]*crhs10 - N[1]*crhs6 + crhs11*(DN(1,0)*crhs7 + DN(1,1)*crhs8);
            rhs[2]=-DN(2,0)*crhs2 - DN(2,1)*crhs4 + N[2]*crhs0 - N[2]*crhs10 - N[2]*crhs6 + crhs11*(DN(2,0)*crhs7 + DN(2,1)*crhs8);

    // KRATOS_WATCH(tau)
}

}
