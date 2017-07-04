//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Christian Rossi
//

#include "heat_equation_Neumann_condition.h"

namespace Kratos
{


template <>
void HeatEquationNeumannCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 2;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(TEMPERATURE).EquationId();


    }
}

/**
 * @see HeatEquationNeumannCondition::EquationIdVector
 */
template <>
void HeatEquationNeumannCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 3;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {

        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(TEMPERATURE).EquationId();
    }
}

/**
 * @see HeatEquationNeumannCondition::GetDofList
 */
template <>
void HeatEquationNeumannCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
                                                ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 2;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {

        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(TEMPERATURE);
    }
}

/**
 * @see HeatEquationNeumannCondition::GetDofList
 */
template <>
void HeatEquationNeumannCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
                                                ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 3;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {

        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(TEMPERATURE);
    }
}

/// Computes the Gauss pt. LHS contribution
/**
* @param lhs_gauss reference to the local LHS matrix
* @param data Gauss pt. data structure
*/
template<unsigned int Tdim, unsigned int TNumNodes>
void HeatEquationNeumannCondition<Tdim,TNumNodes>::ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes,TNumNodes>& lhs_gauss,
const ConditionDataStruct& data)
{
    noalias(lhs_gauss) = ZeroMatrix(TNumNodes,TNumNodes);
}/// Computes the Gauss pt. RHS contribution
/**
* @param rhs_gauss reference to the local RHS vector
* @param data Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void HeatEquationNeumannCondition<TDim,TNumNodes>::ComputeGaussPointRHSContribution(array_1d<double,TNumNodes>& rhs_gauss,
const ConditionDataStruct& data)
{

    // Initialize the local RHS
    noalias(rhs_gauss) = ZeroVector(TNumNodes);

    // Gauss pt. Neumann BC contribution
    this->ComputeRHSNeumannContribution(rhs_gauss, data);
}

/// Computes the condition RHS Neumann BC contribution
/**
* @param rhs_gauss reference to the local RHS vector
* @param data Gauss pt. data structure
*/
template<unsigned int Tdim, unsigned int TNumNodes>
void HeatEquationNeumannCondition<Tdim,TNumNodes>::ComputeRHSNeumannContribution(array_1d<double,TNumNodes>& rhs_gauss,
                                                                              const ConditionDataStruct& data)
{
    const GeometryType& rGeom = this->GetGeometry();

    // Add Neumann BC contribution
    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        const double flux_ext = rGeom[i].FastGetSolutionStepValue(FACE_HEAT_FLUX);
        rhs_gauss[i] -= data.wGauss*data.N[i]*flux_ext;
        std::cout << "Neumann contribution OK\n";
    }
}

/// Computes the 2D condition normal
/**
* @param An reference to condition normal vector
*/
template <>
void HeatEquationNeumannCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;

}

/// Computes the 3D condition normal
/**
* @param An reference to condition normal vector
*/
template <>
void HeatEquationNeumannCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}


template class HeatEquationNeumannCondition<2,2>;
template class HeatEquationNeumannCondition<3,3>;

} // namespace Kratos
